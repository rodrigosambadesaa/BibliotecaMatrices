package comun.resultados;

import comun.MatrizCompleja;

import java.math.BigDecimal;

/** Descomposición singular A=UΣV* con error certificado. */
public record ResultadoSVD<M extends MatrizCompleja>(M u, M sigma, M v, BigDecimal error) {
}
