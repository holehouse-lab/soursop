"""
Tests for soursop.ped (the Protein Ensemble Database client).

The network is replaced by a small fake PED server (``FakePED``) that answers
requests the way the PED API describes (https://proteinensemble.org/api), so
these tests run offline. The ensemble it serves is a real multi-model PDB
built from the bundled GS6 trajectory.

Tests marked ``live`` talk to the real PED server and only run when the
environment variable ``SOURSOP_PED_LIVE=1`` is set.
"""

import copy
import gzip
import io
import json
import os
import tarfile
import urllib.error
import urllib.parse

import mdtraj as md
import numpy as np
import pytest

import soursop
from soursop import ped
from soursop.ped import _http
from soursop.ssexceptions import SSException

test_data_dir = soursop.get_data("test_data")

# Entry record in the shape of the example in PED's API description
ENTRY = {
    "entry_id": "PED00001",
    "version": 1,
    "creation_date": "2020-01-01",
    "description": {
        "title": "Structural ensemble of pSic1 (1-90)",
        "authors": [
            {"name": "Tanja Mittag", "orcid_id": None},
            {"name": "Julie D. Forman-Kay", "orcid_id": None},
        ],
        "publication_identifier": "20399186",
        "entry_cross_reference": [{"db": "disprot", "id": "DP00631"}],
        "experimental_cross_reference": [{"db": "bmrb", "id": "16659"}],
        "ontology_terms": [
            {"name": "NMR", "namespace": "Measurement method", "id": "00120"},
            {"name": "SAXS", "namespace": "Measurement method", "id": "00125"},
        ],
    },
    "construct_chains": [
        {
            "chain_name": "A",
            "fragments": [
                {"description": "tag", "uniprot_acc": None},
                {"description": "Protein SIC1", "uniprot_acc": "P38634"},
            ],
        }
    ],
    "ensembles": [
        {
            "ensemble_id": "e001",
            "models": 5,
            "chains": [
                {
                    "chain_name": "A",
                    "rg_mean": 26.74,
                    "entropy_dssp_mean": 0.40,
                    "relative_asa_mean": 0.69,
                }
            ],
        },
        {
            # the shape PED's live API returns
            "ensemble_id": "e002",
            "models": 10,
            "ensemble_details": {
                "models": 10,
                "only_CA": False,
                "chains": [
                    {"chain_name": "A", "chain_type": "protein", "sequence": "GSMTPS"}
                ],
                "analysis": {
                    "summary_dataframe": [
                        {
                            "ensemble_code": "input_ensemble",
                            "rg_mean": 2.70,
                            "flory_exponent": 0.34,
                        }
                    ]
                },
            },
        },
        {"ensemble_id": "e003", "models": 11, "chains": []},
    ],
}


def _single_ensemble_entry(entry_id):
    entry = copy.deepcopy(ENTRY)
    entry["entry_id"] = entry_id
    entry["ensembles"] = entry["ensembles"][:1]
    return entry


@pytest.fixture(scope="module")
def ensemble_pdb(tmp_path_factory):
    """Bytes of a 5-model PDB file (the GS6 trajectory)."""
    traj = md.load(
        os.path.join(test_data_dir, "gs6_AA.xtc"),
        top=os.path.join(test_data_dir, "gs6_AA.pdb"),
    )
    path = tmp_path_factory.mktemp("ped") / "gs6_models.pdb"
    traj.save_pdb(str(path))
    return path.read_bytes()


class FakePED:
    """Stands in for PED: answers requests by path, and records them."""

    def __init__(self, ensemble_pdb):
        self.requests = []
        self.entries = {
            "PED00001": ENTRY,
            "PED00002": _single_ensemble_entry("PED00002"),
        }
        self.search_results = [
            _single_ensemble_entry(f"PED{i:05d}") for i in range(10, 15)
        ]
        self.pdb_payload = ensemble_pdb
        self.weights_payload = b"model,weight\n1,2\n2,2\n3,2\n4,2\n5,2\n"
        self.assets = {
            "gyration-chains": json.dumps([{"chain": "A", "rg": 26.7}]).encode(),
        }

    def __call__(self, url, timeout):
        parsed = urllib.parse.urlparse(url)
        assert url.startswith(_http.PED_API_URL)
        path = parsed.path[len(urllib.parse.urlparse(_http.PED_API_URL).path) :]
        query = dict(urllib.parse.parse_qsl(parsed.query))
        self.requests.append((path, query))
        parts = [p for p in path.split("/") if p]

        if parts == ["entries"]:
            offset = int(query.get("offset", 0))
            limit = int(query.get("limit", 10))
            page = self.search_results[offset : offset + limit]
            return json.dumps(
                {
                    "count": len(self.search_results),
                    "limit": limit,
                    "offset": offset,
                    "result": page,
                }
            ).encode()
        if len(parts) == 2 and parts[0] == "entries":
            if parts[1] not in self.entries:
                raise _http.PEDError("not found", status=404)
            return json.dumps(self.entries[parts[1]]).encode()
        if len(parts) == 5 and parts[2] == "ensembles":
            asset = parts[4]
            if asset == "ensemble-pdb":
                return self.pdb_payload
            if asset == "weights":
                if self.weights_payload is None:
                    raise _http.PEDError("not found", status=404)
                return self.weights_payload
            if asset in self.assets:
                return self.assets[asset]
        if len(parts) == 5 and parts[4] == "download-all-data":
            return b"\x1f\x8bfake-archive"
        raise _http.PEDError(f"no such resource {path}", status=404)


@pytest.fixture
def fake_ped(monkeypatch, ensemble_pdb):
    server = FakePED(ensemble_pdb)
    monkeypatch.setattr(_http, "_request", server)
    return server


# ------------------------------------------------------------------------------
# identifiers
# ------------------------------------------------------------------------------
class TestIdentifiers:
    def test_parse(self):
        assert ped.parse_identifier("PED00001") == ("PED00001", None)
        assert ped.parse_identifier(" ped00001E001 ") == ("PED00001", "e001")
        assert ped.parse_identifier("PED00001", 2) == ("PED00001", "e002")
        assert ped.parse_identifier("PED00001", "e3") == ("PED00001", "e003")
        assert ped.parse_identifier("PED00001e001", "e001") == ("PED00001", "e001")
        assert ped.full_identifier("PED00001", 1) == "PED00001e001"

    @pytest.mark.parametrize(
        "bad", ["", "PED1", "XYZ00001", "PED00001x001", "PED00001e", 12345]
    )
    def test_malformed_identifiers_raise(self, bad):
        with pytest.raises(SSException):
            ped.parse_identifier(bad)

    def test_conflicting_ensemble_raises(self):
        with pytest.raises(SSException, match="names ensemble e001"):
            ped.parse_identifier("PED00001e001", "e002")

    @pytest.mark.parametrize("bad", [0, -1, "e0", "abc", True])
    def test_bad_ensemble_ids(self, bad):
        with pytest.raises(SSException):
            ped.normalise_ensemble_id(bad)


# ------------------------------------------------------------------------------
# HTTP layer
# ------------------------------------------------------------------------------
class TestHTTP:
    def test_build_url(self):
        url = _http.build_url(
            "/entries", {"free_text": "sic 1", "limit": 5, "x": None, "flag": True}
        )
        assert url.startswith("https://deposition.proteinensemble.org/v1/entries?")
        query = dict(urllib.parse.parse_qsl(urllib.parse.urlparse(url).query))
        assert query == {"free_text": "sic 1", "limit": "5", "flag": "true"}

    def test_unreachable_server_raises_ped_error(self, monkeypatch):
        def fail(*args, **kwargs):
            raise urllib.error.URLError("connection refused")

        monkeypatch.setattr(_http.urllib.request, "urlopen", fail)
        with pytest.raises(ped.PEDError, match="Could not reach PED") as info:
            ped.get_entry("PED00001")
        assert info.value.status is None
        assert isinstance(info.value, SSException)

    def test_http_error_carries_status(self, monkeypatch):
        def not_found(request, timeout):
            raise urllib.error.HTTPError(request.full_url, 404, "Not Found", {}, None)

        monkeypatch.setattr(_http.urllib.request, "urlopen", not_found)
        with pytest.raises(ped.PEDError) as info:
            ped.get_entry("PED99999")
        assert info.value.status == 404

    def test_invalid_json_raises(self, monkeypatch):
        monkeypatch.setattr(
            _http, "_request", lambda url, timeout: b"<html>oops</html>"
        )
        with pytest.raises(ped.PEDError, match="not valid JSON"):
            ped.get_entry("PED00001")


# ------------------------------------------------------------------------------
# entries
# ------------------------------------------------------------------------------
class TestEntries:
    def test_get_entry(self, fake_ped):
        entry = ped.get_entry("PED00001e002")  # ensemble part is ignored
        assert entry.entry_id == "PED00001"
        assert entry.title.startswith("Structural ensemble of pSic1")
        assert entry.version == 1
        assert entry.authors == ("Tanja Mittag", "Julie D. Forman-Kay")
        assert entry.uniprot_accessions == ("P38634",)
        assert entry.protein_names == ("tag", "Protein SIC1")
        assert entry.publication == "20399186"
        assert entry.ontology_terms == ("NMR", "SAXS")
        assert ("disprot", "DP00631") in entry.cross_references
        assert ("bmrb", "16659") in entry.cross_references
        assert entry.ensemble_ids == ("e001", "e002", "e003")
        assert entry.identifiers == ("PED00001e001", "PED00001e002", "PED00001e003")
        e1 = entry.ensemble(1)
        assert e1.n_models == 5 and e1.chains[0].rg_mean == pytest.approx(26.74)
        assert e1.only_ca is None and e1.ped_statistics == {}
        e2 = entry.ensemble("e002")
        assert e2.only_ca is False
        assert (
            e2.chains[0].sequence == "GSMTPS" and e2.chains[0].chain_type == "protein"
        )
        assert e2.ped_statistics["rg_mean"] == pytest.approx(2.70)
        assert e2.ped_statistics["flory_exponent"] == pytest.approx(0.34)
        assert "PED00001" in str(entry)
        assert fake_ped.requests[0][0] == "/entries/PED00001"

    def test_missing_entry(self, fake_ped):
        with pytest.raises(ped.PEDError) as info:
            ped.get_entry("PED99999")
        assert info.value.status == 404

    def test_list_ensembles(self, fake_ped):
        summaries = ped.list_ensembles("PED00001")
        assert [s.identifier for s in summaries] == [
            "PED00001e001",
            "PED00001e002",
            "PED00001e003",
        ]

    def test_unknown_ensemble_in_entry(self, fake_ped):
        with pytest.raises(SSException, match="has no ensemble e009"):
            ped.get_entry("PED00001").ensemble(9)


# ------------------------------------------------------------------------------
# search
# ------------------------------------------------------------------------------
class TestSearch:
    def test_pages_through_all_results(self, fake_ped):
        ids = ped.search("sic1", page_size=2)
        assert ids == ["PED00010", "PED00011", "PED00012", "PED00013", "PED00014"]
        pages = [q for p, q in fake_ped.requests if p == "/entries"]
        assert [int(q["offset"]) for q in pages] == [0, 2, 4]
        assert all(q["free_text"] == "sic1" for q in pages)

    def test_max_results(self, fake_ped):
        assert ped.search("sic1", max_results=3, page_size=2) == [
            "PED00010",
            "PED00011",
            "PED00012",
        ]

    def test_ensemble_identifiers(self, fake_ped):
        ids = ped.search("sic1", ensembles=True, max_results=2)
        assert ids == ["PED00010e001", "PED00011e001"]

    def test_criteria_are_mapped_to_ped_parameters(self, fake_ped):
        ped.search(
            protein_name="Sic1",
            uniprot="P38634",
            term="NMR",
            author="Mittag",
            publication="20399186",
            publication_title="pSic1",
            cross_ref="DP00631",
            max_results=1,
        )
        query = fake_ped.requests[0][1]
        assert query["protein_name"] == "Sic1"
        assert query["uniprot_acc"] == "P38634"
        assert query["term"] == "NMR"
        assert query["data_owner"] == "Mittag"
        assert query["publication_identifier"] == "20399186"
        assert query["publication_html"] == "pSic1"
        assert query["cross_ref"] == "DP00631"
        assert "free_text" not in query

    def test_search_entries_returns_records(self, fake_ped):
        entries = ped.search_entries("sic1", max_results=2)
        assert all(isinstance(e, ped.PEDEntry) for e in entries)
        assert entries[0].identifiers == ("PED00010e001",)

    @pytest.mark.parametrize(
        "kwargs", [{"max_results": 0}, {"page_size": 0}, {"max_results": 1.5}]
    )
    def test_invalid_paging_arguments(self, fake_ped, kwargs):
        with pytest.raises(SSException):
            ped.search("sic1", **kwargs)


# ------------------------------------------------------------------------------
# ensembles
# ------------------------------------------------------------------------------
def _gz(data):
    return gzip.compress(data)


def _tar(data, name="PED00001e001.pdb", compress=True):
    buffer = io.BytesIO()
    mode = "w:gz" if compress else "w"
    with tarfile.open(fileobj=buffer, mode=mode) as tar:
        info = tarfile.TarInfo(name)
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    return buffer.getvalue()


class TestEnsembles:
    @pytest.mark.parametrize(
        "wrap",
        [lambda d: d, _gz, _tar, lambda d: _tar(d, compress=False)],
        ids=["plain", "gzip", "tar.gz", "tar"],
    )
    def test_load_ensemble(self, fake_ped, ensemble_pdb, tmp_path, wrap):
        fake_ped.pdb_payload = wrap(ensemble_pdb)
        traj = ped.load_ensemble("PED00001e001")
        assert traj.n_frames == 5
        reference = tmp_path / "reference.pdb"
        reference.write_bytes(ensemble_pdb)
        np.testing.assert_allclose(traj.traj.xyz, md.load_pdb(str(reference)).xyz)
        rg = traj.proteinTrajectoryList[0].get_radius_of_gyration()
        assert rg.shape == (5,) and np.all(rg > 0)

    def test_entry_and_ensemble_forms_are_equivalent(self, fake_ped):
        a = ped.load_ensemble("PED00001e001")
        b = ped.load_ensemble("PED00001", 1)
        c = ped.get_entry("PED00001").load("e001")
        np.testing.assert_array_equal(a.traj.xyz, b.traj.xyz)
        np.testing.assert_array_equal(a.traj.xyz, c.traj.xyz)

    def test_entry_only_needs_a_single_ensemble(self, fake_ped):
        assert ped.load_ensemble("PED00002").n_frames == 5
        with pytest.raises(SSException, match="contains 3 ensembles"):
            ped.load_ensemble("PED00001")

    def test_cache_dir_downloads_once(self, fake_ped, tmp_path):
        ped.load_ensemble("PED00001e001", cache_dir=tmp_path)
        ped.load_ensemble("PED00001e001", cache_dir=tmp_path)
        downloads = [p for p, q in fake_ped.requests if p.endswith("ensemble-pdb")]
        assert len(downloads) == 1
        assert (tmp_path / "PED00001e001.pdb").exists()

    def test_download_ensemble(self, fake_ped, ensemble_pdb, tmp_path):
        path = ped.download_ensemble("PED00001e001", directory=tmp_path / "out")
        assert path.endswith("PED00001e001.pdb")
        with open(path, "rb") as fh:
            assert fh.read() == ensemble_pdb
        # no temporary files left behind
        assert os.listdir(tmp_path / "out") == ["PED00001e001.pdb"]

    def test_sstrajectory_kwargs_are_passed_through(self, fake_ped):
        traj = ped.load_ensemble("PED00001e001", protein_grouping=[[1, 2, 3, 4, 5, 6]])
        assert traj.proteinTrajectoryList[0].n_residues == 6
        with pytest.raises(SSException, match="sets"):
            ped.load_ensemble("PED00001e001", TRJ=None)

    @pytest.mark.parametrize(
        "payload, match",
        [
            (b"data_PED00001\n_atom_site.id 1\n", "mmCIF"),
            (b"<html>not a structure</html>", "does not look like a PDB"),
            (_tar(b"ATOM  ", name="a.pdb") + b"", None),
        ],
        ids=["mmcif", "garbage", "archive-ok"],
    )
    def test_unreadable_downloads(self, fake_ped, payload, match):
        fake_ped.pdb_payload = payload
        if match is None:
            # a tar with a single (tiny) PDB member is accepted by the extractor
            from soursop.ped.ensemble import _extract_pdb

            assert _extract_pdb(payload, "x") == b"ATOM  "
            return
        with pytest.raises(ped.PEDError, match=match):
            ped.load_ensemble("PED00001e001")

    def test_archive_with_several_pdbs_raises(self):
        from soursop.ped.ensemble import _extract_pdb

        buffer = io.BytesIO()
        with tarfile.open(fileobj=buffer, mode="w") as tar:
            for name in ("a.pdb", "b.pdb"):
                info = tarfile.TarInfo(name)
                info.size = 6
                tar.addfile(info, io.BytesIO(b"ATOM  "))
        with pytest.raises(ped.PEDError, match="expected one PDB file"):
            _extract_pdb(buffer.getvalue(), "x")


class TestWeightsAndAssets:
    def test_text_weights_are_normalised(self, fake_ped):
        w = ped.get_ensemble_weights("PED00001e001")
        np.testing.assert_allclose(w, np.full(5, 0.2))
        raw = ped.get_ensemble_weights("PED00001e001", normalise=False)
        np.testing.assert_allclose(raw, np.full(5, 2.0))

    @pytest.mark.parametrize(
        "payload",
        [
            b"[0.1, 0.2, 0.3, 0.4]",
            b'{"weights": [0.1, 0.2, 0.3, 0.4]}',
            b'[{"model": 1, "weight": 0.1}, {"model": 2, "weight": 0.2},'
            b' {"model": 3, "weight": 0.3}, {"model": 4, "weight": 0.4}]',
            b"# comment\n0.1\n0.2\n0.3\n0.4\n",
        ],
        ids=["json-list", "json-object", "json-records", "text"],
    )
    def test_weight_formats(self, fake_ped, payload):
        fake_ped.weights_payload = payload
        np.testing.assert_allclose(
            ped.get_ensemble_weights("PED00001e001"), [0.1, 0.2, 0.3, 0.4]
        )

    def test_unweighted_ensemble_returns_none(self, fake_ped):
        fake_ped.weights_payload = None
        assert ped.get_ensemble_weights("PED00001e001") is None

    def test_bad_weights_raise(self, fake_ped):
        fake_ped.weights_payload = b"[0.5, -0.1]"
        with pytest.raises(ped.PEDError, match="invalid"):
            ped.get_ensemble_weights("PED00001e001")

    def test_weights_work_with_soursop(self, fake_ped):
        traj = ped.load_ensemble("PED00001e001")
        w = ped.get_ensemble_weights("PED00001e001")
        protein = traj.proteinTrajectoryList[0]
        assert protein.get_radius_of_gyration(weights=w) == pytest.approx(
            protein.get_radius_of_gyration().mean()
        )

    def test_get_ensemble_asset(self, fake_ped):
        data = ped.get_ensemble_asset("PED00001e001", asset="gyration-chains")
        assert data == [{"chain": "A", "rg": 26.7}]
        query = fake_ped.requests[-1][1]
        assert query["response_format"] == "json"
        with pytest.raises(SSException, match="Unknown PED asset"):
            ped.get_ensemble_asset("PED00001e001", asset="not-an-asset")
        with pytest.raises(SSException, match="response_format"):
            ped.get_ensemble_asset(
                "PED00001e001", asset="gyration-chains", response_format="xml"
            )

    def test_download_all_data(self, fake_ped, tmp_path):
        path = ped.download_all_data("PED00001e001", directory=tmp_path)
        assert path.endswith("PED00001e001_all_data.tar.gz")
        assert fake_ped.requests[-1][1]["response_format"] == "json"


# ------------------------------------------------------------------------------
# live tests against the real PED server (opt-in)
# ------------------------------------------------------------------------------
live = pytest.mark.skipif(
    os.environ.get("SOURSOP_PED_LIVE") != "1",
    reason="set SOURSOP_PED_LIVE=1 to run tests against the real PED server",
)


@live
def test_live_entry_and_search():
    entry = ped.get_entry("PED00001")
    assert entry.entry_id == "PED00001" and len(entry.ensembles) >= 1
    assert "PED00001" in ped.search(uniprot="P38634")


@live
def test_live_load_ensemble(tmp_path):
    entry = ped.get_entry("PED00001")
    summary = entry.ensembles[0]
    traj = summary.load(cache_dir=tmp_path)
    if summary.n_models is not None:
        assert traj.n_frames == summary.n_models
    assert traj.proteinTrajectoryList[0].get_radius_of_gyration().shape == (
        traj.n_frames,
    )
