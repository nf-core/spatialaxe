import contextlib
import io
import runpy
import sys
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch


SCRIPT = Path(__file__).resolve().parents[4] / "bin" / "resolvi.py"


class TestResolvi(unittest.TestCase):
    def setUp(self):
        self.scanpy = MagicMock()
        self.pandas = MagicMock()
        self.model_module = MagicMock()
        self.adata = self.scanpy.read_10x_h5.return_value
        self.adata.n_obs = 80
        self.adata.n_vars = 100
        self.modules = {
            "scanpy": self.scanpy,
            "pandas": self.pandas,
            "scviva": MagicMock(),
            "scviva.model": self.model_module,
        }

    def run_script(self, *arguments):
        with patch.dict(sys.modules, self.modules), patch.object(
            sys, "argv", [str(SCRIPT), *arguments]
        ), contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            runpy.run_path(str(SCRIPT), run_name="__main__")

    def test_xenium_input_and_coordinate_alignment(self):
        self.run_script("--input", "xenium", "--output", "result.h5ad")
        self.scanpy.read_10x_h5.assert_called_once_with(Path("xenium/cell_feature_matrix.h5"))
        self.pandas.read_csv.assert_called_once_with(
            Path("xenium/cells.csv.gz"), usecols=["cell_id", "x_centroid", "y_centroid"]
        )
        cells = self.pandas.read_csv.return_value.set_index.return_value
        self.pandas.read_csv.return_value.set_index.assert_called_once_with("cell_id")
        cells.loc.__getitem__.assert_called_once_with(
            (self.adata.obs_names, ["x_centroid", "y_centroid"])
        )
        self.adata.obsm.__setitem__.assert_any_call(
            "spatial", cells.loc.__getitem__.return_value.to_numpy.return_value
        )

    def test_default_training_and_h5ad_output(self):
        self.run_script("--input", "xenium", "--output", "result.h5ad")
        self.scanpy.pp.filter_cells.assert_called_once_with(self.adata, min_counts=5)
        resolvi = self.model_module.ResolVI
        resolvi.setup_anndata.assert_called_once_with(
            self.adata, prepare_data_kwargs={"spatial_rep": "spatial"}
        )
        resolvi.assert_called_once_with(self.adata, semisupervised=False)
        resolvi.return_value.train.assert_called_once_with(max_epochs=100, accelerator="cpu")
        self.adata.obsm.__setitem__.assert_any_call(
            "X_resolVI", resolvi.return_value.get_latent_representation.return_value
        )
        self.adata.write_h5ad.assert_called_once_with("result.h5ad")

    def test_max_epochs_override(self):
        self.run_script("--input", "xenium", "--output", "result.h5ad", "--max-epochs", "1")
        self.model_module.ResolVI.return_value.train.assert_called_once_with(
            max_epochs=1, accelerator="cpu"
        )

    def test_too_few_cells_after_filtering(self):
        def filter_cells(adata, min_counts):
            adata.n_obs = 14

        self.scanpy.pp.filter_cells.side_effect = filter_cells
        with self.assertRaisesRegex(ValueError, "Too few cells"):
            self.run_script("--input", "xenium", "--output", "result.h5ad")
        self.model_module.ResolVI.setup_anndata.assert_not_called()
        self.adata.write_h5ad.assert_not_called()

    def test_required_arguments(self):
        for arguments in [[], ["--input", "xenium"], ["--output", "result.h5ad"]]:
            with self.subTest(arguments=arguments), self.assertRaises(SystemExit) as error:
                self.run_script(*arguments)
            self.assertEqual(error.exception.code, 2)
        self.scanpy.read_10x_h5.assert_not_called()


if __name__ == "__main__":
    unittest.main()