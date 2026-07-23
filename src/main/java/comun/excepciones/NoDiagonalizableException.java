package comun.excepciones;

/** Indica que la matriz no posee una base completa de autovectores. */
public final class NoDiagonalizableException extends OperacionNoDefinidaException {
    public NoDiagonalizableException(String mensaje) { super(mensaje); }
}
