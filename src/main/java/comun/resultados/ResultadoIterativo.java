package comun.resultados;

import comun.MatrizCompleja;

import java.math.BigDecimal;

/** Resultado común de algoritmos iterativos. */
public record ResultadoIterativo<M extends MatrizCompleja>(
        M valor, int iteraciones, BigDecimal residuo, boolean convergio) {
}
