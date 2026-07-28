import tempfile
import unittest
from pathlib import Path

import numpy as np

from ms_deconvolution import parse_input, read_table


class MKSFormatTest(unittest.TestCase):
    def test_quoted_tab_export_with_timestamp_scan_and_trailing_column(self) -> None:
        content = (
            '"Time"\t"Scan"\t"Mass 2"\t"Mass 18"\t"Mass 43"\t"Mass 72"\t\n'
            '"12/4/2025 5:41:48 PM"\t1\t1.0e+02\t3.0e+02\t4.0e+02\t7.2e+01\t\n'
            '"12/4/2025 5:41:53 PM"\t2\t1.1e+02\t3.1e+02\t4.1e+02\t7.3e+01\t\n'
            '"12/4/2025 5:41:58 PM"\t3\t1.2e+02\t3.2e+02\t4.2e+02\t7.4e+01\t\n'
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "mks_export.txt"
            path.write_text(content, encoding="utf-8")
            frame = read_table(path)
            parsed = parse_input(
                frame,
                data_format="auto",
                time_column=None,
                mass_column=None,
                signal_column=None,
                time_unit="auto",
            )

        self.assertEqual(parsed.source_format, "wide MKS scan export")
        self.assertEqual(parsed.source_delimiter, "tab")
        np.testing.assert_allclose(parsed.time_seconds, [0.0, 5.0, 10.0])
        np.testing.assert_allclose(parsed.scan_ids, [1.0, 2.0, 3.0])
        np.testing.assert_allclose(parsed.masses, [2.0, 18.0, 43.0, 72.0])
        self.assertEqual(parsed.original_time[0], "12/4/2025 5:41:48 PM")
        self.assertEqual(parsed.signals.shape, (3, 4))


if __name__ == "__main__":
    unittest.main()
