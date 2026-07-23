package comun.excepciones;

/** Indica que el contexto solicitado no permite certificar el resultado. */
public final class PrecisionInsuficienteException extends MatrizException {
    public PrecisionInsuficienteException(String mensaje) { super(mensaje); }
}
