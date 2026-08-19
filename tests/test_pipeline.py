import tempfile
import unittest
from pathlib import Path

from dreampep.generator import MarkovGenerator
from dreampep.generator import CANONICAL
from dreampep.ingest import ingest_reviewed_csv
from dreampep.io import read_jsonl, write_jsonl
from dreampep.library import build_disulfide_library
from dreampep.pipeline import DesignPipeline
from dreampep.scoring import DevelopabilityScorer, LinearStudent


ROOT = Path(__file__).resolve().parents[1]


class PipelineTest(unittest.TestCase):
    def test_library_design(self):
        rows = build_disulfide_library(24, 10, seed=3)
        self.assertEqual(len(rows), 24)
        self.assertEqual(len({row["sequence"] for row in rows}), 24)
        self.assertTrue(all(row["sequence"].startswith("C") and row["sequence"].endswith("C") for row in rows))

    def test_end_to_end(self):
        records = ingest_reviewed_csv(ROOT / "examples" / "reviewed_patent_examples.csv")
        self.assertEqual(len(records), 4)
        with tempfile.TemporaryDirectory() as directory:
            data = Path(directory) / "data.jsonl"
            write_jsonl(data, records)
            loaded = read_jsonl(data)
            student = LinearStudent()
            student.fit(loaded, epochs=10)
            generator = MarkovGenerator(seed=1)
            generator.fit(loaded)
            pipeline = DesignPipeline(generator, {"activity": student, "dev": DevelopabilityScorer()})
            ranked = pipeline.design("DEMO_TARGET", n=3, min_len=8, max_len=12)
            self.assertTrue(ranked)
            self.assertGreaterEqual(ranked[0].aggregate, ranked[-1].aggregate)
            self.assertTrue(all(set(item.sequence) <= set(CANONICAL) for item in ranked))


if __name__ == "__main__":
    unittest.main()
