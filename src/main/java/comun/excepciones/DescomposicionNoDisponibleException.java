package comun.excepciones;

/** Indica que no se cumplen las hipótesis de una descomposición concreta. */
public final class DescomposicionNoDisponibleException extends OperacionNoDefinidaException {
    public DescomposicionNoDisponibleException(String mensaje) { super(mensaje); }
}
