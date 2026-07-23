package comun.excepciones;

/** Indica que una operación compleja requiere hermiticidad. */
public final class MatrizNoHermiticaException extends OperacionNoDefinidaException {
    public MatrizNoHermiticaException(String mensaje) { super(mensaje); }
}
