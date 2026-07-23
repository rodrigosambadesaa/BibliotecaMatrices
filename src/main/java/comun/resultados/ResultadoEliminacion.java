package comun.resultados;

import comun.MatrizCompleja;

import java.math.BigDecimal;
import java.util.List;

/** Traza inmutable de una eliminación gaussiana. */
public record ResultadoEliminacion<M extends MatrizCompleja>(
        M matriz, List<String> operaciones, List<Integer> permutacionFilas,
        List<int[]> pivotes, int signoIntercambios, List<BigDecimal> factoresEscala,
        BigDecimal tolerancia, int rango, List<Integer> columnasPivote,
        List<Integer> variablesLibres) {
    public ResultadoEliminacion {
        operaciones = List.copyOf(operaciones);
        permutacionFilas = List.copyOf(permutacionFilas);
        pivotes = pivotes.stream().map(int[]::clone).toList();
        factoresEscala = List.copyOf(factoresEscala);
        columnasPivote = List.copyOf(columnasPivote);
        variablesLibres = List.copyOf(variablesLibres);
    }
}
