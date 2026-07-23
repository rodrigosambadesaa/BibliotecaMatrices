package comun.excepciones;

/** Indica que las hipótesis matemáticas de una operación no se cumplen. */
public class OperacionNoDefinidaException extends MatrizException {
    public OperacionNoDefinidaException(String mensaje) { super(mensaje); }
}
