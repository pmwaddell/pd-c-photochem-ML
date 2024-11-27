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
    # configure_logging("test_bigjob")
    #
    # geom_opt_arguments = dict(
    #     functional="B3LYP",
    #     basis_set="def2-TZVP",
    #     newgto='',
    #     RI="RIJCOSX",
    #     dispersion_correction="D3BJ",
    #     solvent="CPCM(Chloroform)",
    #     freq=True
    # )
    #
    # tddft_arguments = dict(
    #     functional="B3LYP",
    #     basis_set="def2-TZVPP",
    #     newgto='',  # LANL2DZ?????????????????????????????????????????????????
    #     RI="RIJCOSX",
    #     dispersion_correction="D3BJ",
    #     solvent="CPCM(Chloroform)"
    # )
    #
    # orca_job_sequence(
    #     path_to_conf_search_xyz_files="data/test_bigjob/conf_search",
    #     destination_path="data/test_bigjob/test_result",
    #     geom_opt_arguments=geom_opt_arguments,
    #     part_2_arguments=tddft_arguments
    # )


    # previous TDDFT job
    # orca_job(
    #     path_to_xyz_file="data/LOT_benchmarks/LOT_benchmark_inputs/8_PdMeCl_def2-TZVP_geom_opt.xyz",
    #     xyz_filename_no_extension="8_PdMeCl_geom_opt",
    #     destination_path="data/LOT_benchmarks/TDDFT/aug-pcseg-2_LANL2DZ",
    #     job_type="TDDFT Calculation",
    #     RI="RIJCOSX",
    #     functional="CAM-B3LYP",
    #     basis_set="aug-pcseg-2",
    #     dispersion_correction="D3BJ",
    #     solvent="CPCM(Chloroform)",
    #     newgto='46 "LANL2DZ"',
    #     grid=""
    # )
    print("doin CAM-B3LYP/def2-TZVPP no LANL2DZ")
    orca_job(
        path_to_xyz_file="data/LOT_benchmarks/optimal_results/8_PdMeCl_geom_opt_def2-TZVP/8_PdMeCl_geom_opt.xyz",
        xyz_filename_no_extension="8_PdMeCl_geom_opt",
        destination_path="data/LOT_benchmarks/TDDFT/def2-TZVPP",
        job_type="TDDFT Calculation",
        RI="RIJCOSX",
        functional="CAM-B3LYP",
        basis_set="def2-TZVPP",
        dispersion_correction="D3BJ",
        solvent="CPCM(Chloroform)",
    )

    print("doin B3LYP/def2-TZVPP no LANL2DZ")
    orca_job(
        path_to_xyz_file="data/LOT_benchmarks/optimal_results/8_PdMeCl_geom_opt_def2-TZVP/8_PdMeCl_geom_opt.xyz",
        xyz_filename_no_extension="8_PdMeCl_geom_opt",
        destination_path="data/LOT_benchmarks/TDDFT/B3LYP_def2-TZVPP",
        job_type="TDDFT Calculation",
        RI="RIJCOSX",
        functional="B3LYP",
        basis_set="def2-TZVPP",
        dispersion_correction="D3BJ",
        solvent="CPCM(Chloroform)"
    )

    print("doin B3LYP/def2-TZVPP with LANL2DZ")
    orca_job(
        path_to_xyz_file="data/LOT_benchmarks/optimal_results/8_PdMeCl_geom_opt_def2-TZVP/8_PdMeCl_geom_opt.xyz",
        xyz_filename_no_extension="8_PdMeCl_geom_opt",
        destination_path="data/LOT_benchmarks/TDDFT/B3LYP_def2-TZVPP",
        job_type="TDDFT Calculation",
        RI="RIJCOSX",
        functional="B3LYP",
        basis_set="def2-TZVPP",
        newgto='46 "LANL2DZ"',
        dispersion_correction="D3BJ",
        solvent="CPCM(Chloroform)"
    )
