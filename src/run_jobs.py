import logging

# from networkx.algorithms.centrality import dispersion

logger = logging.getLogger(__name__)

from orca_jobs import orca_job_sequence, orca_job


def configure_logging(log_filename: str="src/logs/log.log") -> None:
    """Configure the logger for our purposes."""
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(message)s',
        datefmt='%Y/%m/%d %H:%M:%S',
        handlers=[logging.FileHandler(f"src/logs/{log_filename}.log", mode='w'),
                  stream_handler]
    )


if __name__ == "__main__":
    configure_logging("round_3_pt1")
    geom_opt_arguments = dict(
        functional="B3LYP",
        basis_set="def2-TZVP",
        newgto='',
        RI="RIJCOSX",
        dispersion_correction="D3BJ",
        solvent="CPCM(Chloroform)",
        grid="",
        freq=True,
        NMR=False,
        orbitals=False
    )

    pt_two_arguments = dict(
        functional="B3LYP",
        basis_set="def2-TZVPP",
        newgto='',
        RI="RIJCOSX",
        dispersion_correction="D3BJ",
        solvent="CPCM(Chloroform)",
        grid="",
        freq=False,
        NMR=False,
        orbitals=False
    )

    orca_job_sequence(
        path_to_conf_search_xyz_files="src/data/all_bigjob/round_3/input_files",
        destination_path="src/data/all_bigjob/round_3/round_3_result",
        geom_opt_arguments=geom_opt_arguments,
        part_2_arguments=pt_two_arguments,
        geom_opt=True,
        single_pt=True,
        tddft=True
    )
