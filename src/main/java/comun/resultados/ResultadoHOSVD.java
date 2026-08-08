package comun.resultados;

import comun.Matriz3DCompleja;
import comun.MatrizCompleja;

import java.math.BigDecimal;
import java.util.Objects;

/** Resultado inmutable de una descomposición Tucker/HOSVD de orden tres. */
public record ResultadoHOSVD<T extends Matriz3DCompleja, P extends MatrizCompleja>(
        T nucleo, P factorCapas, P factorFilas, P factorColumnas, BigDecimal error) {
    public ResultadoHOSVD {
        Objects.requireNonNull(nucleo, "El núcleo no puede ser null");
        Objects.requireNonNull(factorCapas, "El factor de capas no puede ser null");
        Objects.requireNonNull(factorFilas, "El factor de filas no puede ser null");
        Objects.requireNonNull(factorColumnas, "El factor de columnas no puede ser null");
        Objects.requireNonNull(error, "El error no puede ser null");
        if (error.signum() < 0) throw new IllegalArgumentException("El error no puede ser negativo");
    }
}
