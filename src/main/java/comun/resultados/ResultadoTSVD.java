package comun.resultados;

import comun.Matriz3DCompleja;

import java.math.BigDecimal;
import java.util.Objects;

/** Descomposición singular asociada al t-producto: A = U * S * Vᴴ. */
public record ResultadoTSVD<T extends Matriz3DCompleja>(T u, T sigma, T v, BigDecimal error) {
    public ResultadoTSVD {
        Objects.requireNonNull(u, "U no puede ser null");
        Objects.requireNonNull(sigma, "Sigma no puede ser null");
        Objects.requireNonNull(v, "V no puede ser null");
        Objects.requireNonNull(error, "El error no puede ser null");
        if (error.signum() < 0) throw new IllegalArgumentException("El error no puede ser negativo");
    }
}
