package comun.excepciones;

/** Indica dimensiones negativas, irregulares o incompatibles. */
public final class DimensionMatrizException extends MatrizException {
    public DimensionMatrizException(String mensaje) { super(mensaje); }
}
