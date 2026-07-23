package comun.excepciones;

/** Indica que se requiere una matriz invertible. */
public final class MatrizSingularException extends MatrizException {
    public MatrizSingularException(String mensaje) { super(mensaje); }
}
