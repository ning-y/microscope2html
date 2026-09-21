import os
import tempfile
import unittest
from pathlib import Path

from PIL import Image

from microscope2html.evos import discover_scan
from microscope2html.pipeline import (
    _translated_registered_config,
    create_tiled_html_channels,
    generate_macro,
)


class MultiChannelTests(unittest.TestCase):
    def test_discovery_keeps_all_channels_by_default(self):
        with tempfile.TemporaryDirectory() as directory:
            for channel in (0, 1):
                for field in (0, 1):
                    Path(directory, (
                        f"scan_Slide_R_p00_0_A01f{field:02d}d{channel}.TIF"
                    )).touch()

            units, problems = discover_scan(directory)

            self.assertEqual(problems, [])
            self.assertEqual(len(units), 1)
            self.assertEqual(sorted(units[0].raw_by_channel), [0, 1])
            self.assertEqual(len(units[0].raw_by_channel[0]), 2)
            self.assertEqual(len(units[0].raw_by_channel[1]), 2)

    def test_explicit_channel_still_filters_discovery(self):
        with tempfile.TemporaryDirectory() as directory:
            for channel in (0, 1):
                Path(directory, (
                    f"scan_Slide_R_p00_0_A01f00d{channel}.TIF"
                )).touch()

            units, problems = discover_scan(directory, channel=1)

            self.assertEqual(problems, [])
            self.assertEqual(list(units[0].raw_by_channel), [1])

    def test_registered_config_is_retargeted_to_sibling_channel(self):
        reference = [("/data/field00d0.TIF", 1, 0, 0, 1, 1)]
        sibling = [("/data/field00d1.TIF", 1, 0, 0, 1, 1)]
        config = "field00d0.TIF; ; (12.5, 9.5)\n"

        translated = _translated_registered_config(
            config, reference, sibling, "/work")

        self.assertIn("../data/field00d1.TIF", translated)
        self.assertNotIn("field00d0.TIF", translated)

    def test_reused_registration_disables_overlap_computation(self):
        with tempfile.TemporaryDirectory() as directory:
            macro = os.path.join(directory, "stitch.ijm")
            generate_macro(macro, "registered.txt", "out.tif", directory,
                           False, compute_overlap=False)

            text = Path(macro).read_text()

            self.assertNotIn("compute_overlap", text)
            self.assertIn("computation_parameters", text)

    def test_viewer_uses_named_additive_layers_and_fading_controls(self):
        with tempfile.TemporaryDirectory() as directory:
            dapi = os.path.join(directory, "dapi.png")
            gfp = os.path.join(directory, "gfp.png")
            output = os.path.join(directory, "viewer.html")
            Image.new("RGB", (2, 2), (0, 0, 64)).save(dapi)
            Image.new("RGB", (2, 2), (0, 64, 0)).save(gfp)

            create_tiled_html_channels({0: dapi, 1: gfp}, output, 1.0)
            html = Path(output).read_text()

            self.assertIn("compositeOperation: 'lighter'", html)
            self.assertIn("success: event => { layerItems[key] = event.item; }", html)
            self.assertIn("order.forEach(key => layerItems[key].setOpacity", html)
            self.assertIn("top: 10px; right: 10px", html)
            self.assertIn("hideControlsTimer", html)
            self.assertIn('data-channel="merge"', html)


if __name__ == "__main__":
    unittest.main()
