import org.junit.jupiter.api.Test;
import otros.MatrizBooleana;
import otros.MatrizModular;

import static org.junit.jupiter.api.Assertions.*;

class OtrosCuerposTest {
    @Test void booleanas() {
        MatrizBooleana a = new MatrizBooleana(new boolean[][]{{false, true}, {false, false}});
        assertTrue(a.clausuraTransitiva().get(0, 1));
        assertEquals(a, a.suma(a));
        assertFalse(a.producto(a).get(0, 1));
        boolean[][] copia = a.getMatriz(); copia[0][1] = false;
        assertTrue(a.get(0, 1));
    }
    @Test void modulares() {
        MatrizModular a = new MatrizModular(new int[][]{{1, 2}, {3, 4}}, 5);
        assertEquals(3, a.determinante());
        assertEquals(2, a.rango());
        assertEquals(new MatrizModular(new int[][]{{1, 0}, {0, 1}}, 5), a.producto(a.inversa()));
        assertEquals(new MatrizModular(new int[][]{{2, 4}, {1, 3}}, 5), a.suma(a));
        assertThrows(RuntimeException.class, () -> new MatrizModular(new int[][]{{1}}, 4));
    }
}
