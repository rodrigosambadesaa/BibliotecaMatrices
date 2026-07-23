package comun.resultados;

import comun.MatrizCompleja;

import java.math.BigDecimal;

/** Factores Q y R con error de reconstrucción y ortogonalidad. */
public record ResultadoQR<M extends MatrizCompleja>(M q, M r, BigDecimal error,
                                                    BigDecimal errorOrtogonalidad) {
}
