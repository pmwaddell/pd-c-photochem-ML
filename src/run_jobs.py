import logging
logger = logging.getLogger(__name__)

from orca_jobs import orca_job_sequence


def configure_logging(log_filename: str="logs/log.log") -> None:
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(message)s',
        datefmt='%H:%M:%S',
        handlers=[logging.FileHandler("logs/log.log", mode='w'),
                  stream_handler]
    )


if __name__ == "__main__":
    configure_logging()
    # Apparently, B3LYP is fine for geometry optimizations here.
    geom_opt_arguments = dict(
        functional="B3LYP",
        basis_set="6-311++G(d,p)",
        newgto='46 "LANL2DZ"',
        RI="RIJCOSX",
        dispersion_correction="D3BJ",
        solvent="CPCM(Chloroform)"
    )

    single_pt_arguments = dict(
        functional="CAM-B3LYP",
        basis_set="6-311++G(d,p)",
        newgto='46 "LANL2DZ"',
        RI="RIJCOSX",
        dispersion_correction="D3BJ",
        solvent="CPCM(Chloroform)",
        NMR=True,
        freq=True
    )

    orca_job_sequence(
        path_to_conf_search_xyz_files='data/dft_test',
        destination_path='data/batch_1_CAM_B3LYP',
        geom_opt_arguments=geom_opt_arguments,
        single_pt_arguments=single_pt_arguments
    )
