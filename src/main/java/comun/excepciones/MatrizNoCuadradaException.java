package comun.excepciones;

/** Indica que una operación requiere una matriz cuadrada. */
public final class MatrizNoCuadradaException extends MatrizException {
    public MatrizNoCuadradaException(String mensaje) { super(mensaje); }
}
