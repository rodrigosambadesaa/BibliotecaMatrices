package comun.excepciones;

/** Indica un índice de fila o columna fuera de rango. */
public final class IndiceMatrizException extends MatrizException {
    public IndiceMatrizException(String mensaje) { super(mensaje); }
}
