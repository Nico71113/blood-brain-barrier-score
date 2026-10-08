import base64
import csv
import io
import math
import unittest
from rdkit import Chem
from bbb.chemistry import molecule_data
from bbb.pipeline import analyze, batch_analyze, read_file, results_csv, uploaded_records
from bbb.scoring import score_descriptors


class ScoringTests(unittest.TestCase):
    def setUp(self):
        # Public author-calculator example, not a research-table fixture.
        self.d = dict(AroR=3, HA=28, MW=368.52, HBA=2, HBD=0, TPSA=6.48, pKa=8.1)

    def test_author_calculator_example(self):
        # Author's Excel example; independent cached result from C11.
        d = dict(AroR=3, HA=28, MW=368.52, HBA=2, HBD=0, TPSA=6.48, pKa=8.1)
        result = score_descriptors(d)
        self.assertAlmostEqual(result["bbb_score"], 5.040301056516096, places=12)
        self.assertAlmostEqual(result["bbb_score"], sum(c["contribution"] for c in result["contributions"]))

    def test_exact_branch_boundaries(self):
        for field, boundary in [("HA", 5), ("TPSA", 0), ("pKa", 3)]:
            result=score_descriptors({**self.d, field: boundary})
            self.assertEqual(next(c["p"] for c in result["contributions"] if c["name"]==field),0)
        for field, boundary in [("HA", 45), ("TPSA", 120), ("pKa", 11)]:
            result=score_descriptors({**self.d, field: boundary})
            self.assertNotEqual(next(c["p"] for c in result["contributions"] if c["name"]==field),0)
            above=score_descriptors({**self.d, field: boundary+1})
            self.assertEqual(next(c["p"] for c in above["contributions"] if c["name"]==field),0)
        d={**self.d, "MW":400, "HBA":1, "HBD":0}
        self.assertEqual(next(c["p"] for c in score_descriptors(d)["contributions"] if c["name"]=="MWHBN"),0)
        d["HBA"]=9
        self.assertGreater(next(c["p"] for c in score_descriptors(d)["contributions"] if c["name"]=="MWHBN"),0)
        d["HBA"]=10
        self.assertEqual(next(c["p"] for c in score_descriptors(d)["contributions"] if c["name"]=="MWHBN"),0)

    def test_missing_zero_and_multivalue_pka(self):
        self.assertEqual(analyze({**self.d,"pKa":None})["score_status"],"missing_input")
        self.assertEqual(analyze({**self.d,"pKa":0})["score_status"],"ok")
        self.assertEqual(analyze({**self.d,"pKa":"2.81\n6"})["score_status"],"invalid_input")

    def test_nonfinite_and_invalid_counts(self):
        for value in [math.nan, math.inf, True, "NaN"]:
            with self.assertRaises(ValueError): score_descriptors({**self.d,"MW":value})
        for value in [-1, 2.5]:
            with self.assertRaises(ValueError): score_descriptors({**self.d,"HA":value})


class StructureTests(unittest.TestCase):
    def test_same_formula_different_graphs(self):
        ethanol=molecule_data("CCO")
        ether=molecule_data("COC")
        self.assertEqual(ethanol["formula"],"C2H6O")
        self.assertEqual(ethanol["formula"],ether["formula"])
        self.assertNotEqual(ethanol["canonical_smiles"],ether["canonical_smiles"])
        self.assertNotEqual(ethanol["svg"],ether["svg"])

    def test_2d_coordinates_and_stereochemical_roundtrip(self):
        data=molecule_data("N[C@@H](C)C(=O)O")
        mol=Chem.MolFromMolBlock(data["molblock"])
        self.assertFalse(mol.GetConformer().Is3D())
        self.assertTrue(all(abs(p[2])<1e-8 for p in mol.GetConformer().GetPositions()))
        self.assertEqual(Chem.MolToSmiles(mol,isomericSmiles=True),data["canonical_smiles"])

    def test_invalid_structures(self):
        for smiles in ["CC(", "CO(C)C", "CCO.[Na+]", "*CC", ""]:
            with self.assertRaises(ValueError): molecule_data(smiles)

    def test_unknown_smiles_not_demo_lookup(self):
        data=molecule_data("CC1CCC(O)CC1")
        self.assertEqual(data["formula"],"C7H14O")
        self.assertIn("<svg",data["svg"])

    def test_pka_never_implicitly_guessed(self):
        no_pka=analyze({"smiles":"CCO"})
        self.assertEqual(no_pka["structure_status"],"ok")
        self.assertIsNone(no_pka["bbb_score"])
        neutral=analyze({"smiles":"CCO", "neutral":True})
        self.assertEqual(neutral["descriptors"]["pKa"],8.81)
        self.assertEqual(neutral["pka_source"],"user_confirmed_neutral_convention")
        self.assertEqual(analyze({"smiles":"CCO", "neutral":True,"pKa":4})["score_status"],"invalid_input")


class ImportTests(unittest.TestCase):
    def test_original_headers_and_descriptor_only_scoring(self):
        row={"compund ID":"public_calculator_example", "Number of Aromatic Rings (Aro_R)":3,
             "Number of Heavy Atoms (HA)":28,"Molecular Weight (MW)":368.52,
             "Number of Hydrogen Bond Acceptor (HBA)":2,"Number of Hydrogen Bond Donor (HBD)":0,
             "Topological Polar Surface Area(TPSA)":6.48,"pKa":8.1,"BBB SCORE":"5.0403"}
        result=analyze(row)
        self.assertEqual(round(result["bbb_score"],4),5.0403)
        self.assertEqual(result["bbb_score_input"],"5.0403")
        self.assertEqual(result["structure_status"],"missing_structure")
        self.assertEqual(analyze({**row,"compund ID":0})["compound_id"],"0")

    def test_repeated_ids_preserve_order_and_row_errors(self):
        results=batch_analyze([{"compound_id":"same","smiles":"CCO","pKa":8},
                               {"compound_id":"same","smiles":"COC","pKa":"3\n4"}])["results"]
        self.assertEqual([r["row_id"] for r in results],[1,2])
        self.assertEqual([r["compound_id"] for r in results],["same","same"])
        self.assertEqual(results[1]["score_status"],"invalid_input")
        self.assertIsNone(results[0]["structure"]["svg"])

    def test_csv_quotes_and_duplicate_headers(self):
        data=b'compound_id,SMILES,pKa\n"name, quoted",CCO,8.1\n'
        rows=uploaded_records({"filename":"test.csv","content_base64":base64.b64encode(data).decode()})
        self.assertEqual(rows[0]["compound_id"],"name, quoted")
        with self.assertRaises(ValueError): read_file("x.csv",b"SMILES,SMILES\nCCO,COC\n")

    def test_excel_formula_rejected_not_using_cached_result(self):
        # Synthetic OOXML fixture: no spreadsheet artifact is authored here.
        import zipfile
        stream=io.BytesIO()
        files={
            "[Content_Types].xml": '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/><Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/></Types>',
            "_rels/.rels": '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>',
            "xl/workbook.xml": '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="Sheet1" sheetId="1" r:id="rId1"/></sheets></workbook>',
            "xl/_rels/workbook.xml.rels": '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/></Relationships>',
            "xl/worksheets/sheet1.xml": '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData><row r="1"><c r="A1" t="inlineStr"><is><t>pKa</t></is></c></row><row r="2"><c r="A2"><f>8+1</f><v>9</v></c></row></sheetData></worksheet>',
        }
        with zipfile.ZipFile(stream,"w") as z:
            for name,data in files.items(): z.writestr(name,data)
        with self.assertRaisesRegex(ValueError,"formulas"): read_file("x.xlsx",stream.getvalue())

    def test_csv_export_injection_and_missing_score(self):
        rows=[analyze({"compound_id":"=1+1","smiles":"CCO"})]
        result=list(csv.DictReader(io.StringIO(results_csv(rows))))[0]
        self.assertEqual(result["compound_id"],"'=1+1")
        self.assertEqual(result["bbb_score"],"")
        self.assertEqual(result["score_status"],"missing_input")


if __name__ == "__main__":
    unittest.main()
