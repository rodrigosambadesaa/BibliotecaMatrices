package comun.excepciones;

/** Indica que un algoritmo iterativo agotó el máximo de iteraciones. */
public final class ConvergenciaException extends MatrizException {
    public ConvergenciaException(String mensaje) { super(mensaje); }
}
