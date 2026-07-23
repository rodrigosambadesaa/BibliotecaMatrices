package comun.excepciones;

/** Indica que una operación espectral requiere normalidad. */
public final class MatrizNoNormalException extends OperacionNoDefinidaException {
    public MatrizNoNormalException(String mensaje) { super(mensaje); }
}
