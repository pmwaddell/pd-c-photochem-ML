from src.orca_jobs import orca_batch_job


if __name__ == "__main__":
    print("\n\nNi ethylene complex single point calcs: ")
    orca_batch_job(path_to_xyz_files="data/diimine_data/geom_opt_NiE",
                   destination_path="data/diimine_data/single_pt_NiE",
                   functional="B3LYP",
                   basis_set="def2-TZVP",
                   RI='RIJCOSX',
                   dispersion_correction="D3BJ",
                   charge=1,
                   NMR=True,
                   freq=True,
                   job_type="Single Point Calculation")
