package comun.excepciones;

/** Base de las excepciones de dominio de la biblioteca matricial. */
public class MatrizException extends RuntimeException {
    public MatrizException(String mensaje) {
        super(mensaje);
    }
}
