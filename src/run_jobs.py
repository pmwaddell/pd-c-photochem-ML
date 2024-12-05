import logging

from networkx.algorithms.centrality import dispersion

logger = logging.getLogger(__name__)

from orca_jobs import orca_job_sequence, orca_job


def configure_logging(log_filename: str="logs/log.log") -> None:
    """Configure the logger for our purposes."""
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(message)s',
        datefmt='%H:%M:%S',
        handlers=[logging.FileHandler(f"logs/{log_filename}.log", mode='w'),
                  stream_handler]
    )


if __name__ == "__main__":
    configure_logging("bigjob_split1_2")

    geom_opt_arguments = dict(
        functional="B3LYP",
        basis_set="def2-TZVP",
        newgto='',
        RI="RIJCOSX",
        dispersion_correction="D3BJ",
        solvent="CPCM(Chloroform)",
        grid="",
        freq=True,
        NMR=False
    )

    tddft_arguments = dict(
        functional="B3LYP",
        basis_set="def2-TZVPP",
        newgto='',
        RI="RIJCOSX",
        dispersion_correction="D3BJ",
        solvent="CPCM(Chloroform)",
        grid="",
        freq=False,
        NMR=False
    )

    orca_job_sequence(
        path_to_conf_search_xyz_files="data/bigjob/conf_search_splits/split1",
        destination_path="data/bigjob/split1_result",
        geom_opt_arguments=geom_opt_arguments,
        part_2_arguments=tddft_arguments,
        tddft=True
    )
