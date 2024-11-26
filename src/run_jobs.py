import logging
logger = logging.getLogger(__name__)

from orca_jobs import orca_job, make_uv_vis_plot


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
    configure_logging("lot_benchmarks")

    # TDDFT time
    # orca_job(
    #     path_to_xyz_file="data/LOT_benchmarks/LOT_benchmark_inputs/8_PdMeCl_def2-TZVP_geom_opt.xyz",
    #     xyz_filename_no_extension="8_PdMeCl_geom_opt",
    #     destination_path="data/LOT_benchmarks/TDDFT/pcseg-2_LANL2DZ",
    #     job_type="TDDFT Calculation",
    #     RI="autoaux RIJCOSX",
    #     functional="CAM-B3LYP",
    #     basis_set="pcseg-2",
    #     dispersion_correction="D3BJ",
    #     solvent="CPCM(Chloroform)",
    #     newgto='46 "LANL2DZ"',
    #     grid=""
    # )

    orca_job(
        path_to_xyz_file="data/LOT_benchmarks/LOT_benchmark_inputs/8_PdMeCl_def2-TZVP_geom_opt.xyz",
        xyz_filename_no_extension="8_PdMeCl_geom_opt",
        destination_path="data/LOT_benchmarks/TDDFT/aug-pcseg-2_LANL2DZ",
        job_type="TDDFT Calculation",
        RI="RIJCOSX",
        functional="CAM-B3LYP",
        basis_set="aug-pcseg-2",
        dispersion_correction="D3BJ",
        solvent="CPCM(Chloroform)",
        newgto='46 "LANL2DZ"',
        grid=""
    )

    # orca_job(
    #     path_to_xyz_file="data/LOT_benchmarks/LOT_benchmark_inputs/8_PdMeCl_def2-TZVP_geom_opt.xyz",
    #     xyz_filename_no_extension="8_PdMeCl_geom_opt",
    #     destination_path="data/LOT_benchmarks/TDDFT/6311",
    #     job_type="TDDFT Calculation",
    #     RI="autoaux RIJCOSX",
    #     functional="CAM-B3LYP",
    #     basis_set="6-311G++G(d,p)",
    #     dispersion_correction="D3BJ",
    #     solvent="CPCM(Chloroform)",
    #     newgto='',
    #     grid=""
    # )
