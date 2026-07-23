package comun.resultados;

import ch.obermuhlner.math.big.BigComplex;
import comun.MatrizCompleja;

import java.math.BigDecimal;

/** Autopar numérico con convergencia y residuo. */
public record ResultadoAutovalor<M extends MatrizCompleja>(
        BigComplex autovalor, M autovector, int iteraciones,
        BigDecimal residuo, boolean convergio) {
}
