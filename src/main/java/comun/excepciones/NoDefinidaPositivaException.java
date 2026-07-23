package comun.excepciones;

/** Indica que una factorización requiere definitud positiva. */
public final class NoDefinidaPositivaException extends OperacionNoDefinidaException {
    public NoDefinidaPositivaException(String mensaje) { super(mensaje); }
}
