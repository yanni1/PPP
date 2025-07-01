# Animation update function
import Calc_CPP # pyright: ignore[reportMissingImports]
import params_module # pyright: ignore[reportMissingImports]
import Genetic_algo # pyright: ignore[reportMissingImports]

def get_best_params(nt):
    best_tau, best_rho = Genetic_algo.Genetic_algo(nt)
    tau = float(best_tau)
    rho = float(best_rho)
    return best_tau, best_rho