package comun.excepciones;

/** Indica una división entre cero exacto o numérico según la tolerancia. */
public final class DivisionPorCeroException extends MatrizException {
    public DivisionPorCeroException(String mensaje) { super(mensaje); }
}
