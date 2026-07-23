import org.junit.jupiter.api.Test;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Ejecuta la auditoría reflectiva del principal: cada constructor y método
 * público se invoca individualmente y cualquier ausencia hace fallar la prueba.
 */
class CoberturaApiTest {
    @Test void todosLosConstructoresYMetodosPublicosSeDemuestran() {
        PrintStream anterior = System.out;
        ByteArrayOutputStream salida = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(salida, true, StandardCharsets.UTF_8));
            BibliotecaMatricesApp.main(new String[0]);
        } finally {
            System.setOut(anterior);
        }
        String texto = salida.toString(StandardCharsets.UTF_8);
        assertTrue(texto.contains("fallidos=0, no probados=0"), texto);
        assertTrue(texto.contains("array.Matriz#suma"), "No se auditó la API array");
        assertTrue(texto.contains("lista.Matriz#suma"), "No se auditó la API lista");
    }
}
