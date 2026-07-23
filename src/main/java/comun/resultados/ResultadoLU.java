package comun.resultados;

import comun.MatrizCompleja;

import java.math.BigDecimal;

/** Factores P, L y U con error de reconstrucción. */
public record ResultadoLU<M extends MatrizCompleja>(M p, M l, M u, int signo, BigDecimal error) {
}
