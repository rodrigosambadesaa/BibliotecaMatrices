package comun.resultados;

import ch.obermuhlner.math.big.BigComplex;
import comun.MatrizCompleja;

import java.math.BigDecimal;
import java.util.List;

/** Diagonalización A=PDP^-1 y sus residuos. */
public record ResultadoDiagonalizacion<M extends MatrizCompleja>(
        M p, M d, M inversaP, List<BigComplex> autovalores,
        M autovectores, BigDecimal residuo, boolean verificada) {
    public ResultadoDiagonalizacion {
        autovalores = List.copyOf(autovalores);
    }
}
