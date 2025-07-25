# Animation update function
import Calc_CPP # pyright: ignore[reportMissingImports]
import params_module # pyright: ignore[reportMissingImports]
import Genetic_algo # pyright: ignore[reportMissingImports]

def get_best_params(p):
    best_tau, best_rho = Genetic_algo.Genetic_algo(p)
    return int(best_tau), float(best_rho)