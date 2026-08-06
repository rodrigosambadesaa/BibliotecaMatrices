import ch.obermuhlner.math.big.BigComplex;
import comun.Complejos;
import comun.excepciones.DimensionMatrizException;
import comun.excepciones.DivisionPorCeroException;
import comun.excepciones.IndiceMatrizException;
import org.junit.jupiter.api.Test;

import java.lang.reflect.Field;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

class Matriz3DTest {
    private static final BigDecimal TOL = new BigDecimal("1E-20");
    private static final double[][][] DATOS = {
            {{1, 2, 3}, {4, 5, 6}},
            {{7, 8, 9}, {10, 11, 12}}
    };

    @Test void construccionAccesoInmutabilidadYAlmacenamientos() throws Exception {
        double[][][] origen = {
                {{1, 2}, {3, 4}},
                {{5, 6}, {7, 8}}
        };
        array.Matriz3D a = new array.Matriz3D(origen);
        origen[0][0][0] = 99;
        assertEquals(BigComplex.ONE, a.get(0, 0, 0));
        BigComplex[][][] copia = a.getMatriz();
        copia[0][0][0] = BigComplex.ZERO;
        assertEquals(BigComplex.ONE, a.get(0, 0, 0));
        assertEquals("2x2x2", a.dimensiones());
        assertEquals(8, a.numeroElementos());
        assertThrows(IndiceMatrizException.class, () -> a.get(2, 0, 0));
        assertThrows(DimensionMatrizException.class, () -> new array.Matriz3D(new double[][][]{{{1}}, {{1, 2}}}));

        lista.Matriz3D l = new lista.Matriz3D(a);
        assertEquals(a, l);
        assertEquals(a, l.comoMatriz3DArray());
        assertEquals(l, a.comoMatriz3DLista());
        assertNotSame(l.comoListaProfunda(), l.comoListaProfunda());
        assertInstanceOf(lista.Matriz3D.class, l.suma(l));

        Field campoArray = array.Matriz3D.class.getDeclaredField("matriz");
        Field campoLista = lista.Matriz3D.class.getDeclaredField("matriz");
        assertEquals(BigComplex[][][].class, campoArray.getType());
        assertEquals(ArrayList.class, campoLista.getType());

        array.Matriz3D unidad = new array.Matriz3D(new double[][][]{{{1}}});
        assertEquals(unidad, new array.Matriz3D(new Double[][][]{{{1.0}}}));
        assertEquals(unidad, new array.Matriz3D(new BigDecimal[][][]{{{BigDecimal.ONE}}}));
        assertEquals(unidad, new array.Matriz3D(new String[][][]{{{"1"}}}));
        assertEquals(unidad, new array.Matriz3D(new BigComplex[][][]{{{BigComplex.ONE}}}));
        assertEquals(unidad, lista.Matriz3D.desdeNumeros(List.of(List.of(List.of(1)))));
        assertEquals(unidad, lista.Matriz3D.desdeBigComplex(List.of(List.of(List.of(BigComplex.ONE)))));
        assertEquals(MathContext.DECIMAL64, new array.Matriz3D(1, 1, 1, MathContext.DECIMAL64).mathContext());
    }

    @Test void cortesEdicionExtraccionYCapas() {
        array.Matriz3D a = new array.Matriz3D(DATOS);
        assertEquals(new array.Matriz(new double[][]{{1, 2, 3}, {4, 5, 6}}), a.capa(0));
        assertEquals(new array.Matriz(new double[][]{{1, 2, 3}, {7, 8, 9}}), a.cortePorFila(0));
        assertEquals(new array.Matriz(new double[][]{{3, 6}, {9, 12}}), a.cortePorColumna(2));
        assertEquals(99, a.conElemento(1, 1, 2, Complejos.de(99)).get(1, 1, 2).re.intValueExact());
        assertEquals(a, a.conCapa(0, a.capa(0)));
        assertEquals("1x1x2", a.submatriz(1, 2, 1, 2, 1, 3).dimensiones());
        assertEquals(a.capa(1), a.extraerCapas(1).capa(0));
        assertEquals(a, a.eliminarCapa(0).insertarCapa(0, a.capa(0)));
        assertThrows(DimensionMatrizException.class, () -> a.eliminarCapa(0).eliminarCapa(0));
    }

    @Test void aritmeticaNormasYReducciones() {
        array.Matriz3D a = new array.Matriz3D(new double[][][]{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}});
        array.Matriz3D cero = array.Matriz3D.cero(2, 2, 2);
        assertEquals(a, a.suma(cero));
        assertTrue(a.resta(a).esNula());
        assertEquals(a.productoEscalar(2), a.suma(a));
        assertEquals(a, a.productoHadamard(array.Matriz3D.rellena(2, 2, 2, BigComplex.ONE)));
        assertEquals(Complejos.de(36), a.sumaElementos());
        assertEquals(Complejos.de("4.5"), a.promedio());
        assertEquals(new array.Matriz(new double[][]{{6, 8}, {10, 12}}), a.sumaPorCapas());
        assertEquals(new array.Matriz(new double[][]{{3, 4}, {5, 6}}), a.promedioPorCapas());
        assertEquals(Complejos.de(204), a.productoInternoFrobenius(a));
        assertEquals(a.productoInternoFrobenius(a), a.productoFrobenius(a));
        assertEquals(0, a.normaFrobenius().compareTo(new BigDecimal("204").sqrt(a.mathContext())));
        assertEquals(0, a.normaMaximo().compareTo(BigDecimal.valueOf(8)));
        assertThrows(DivisionPorCeroException.class, () -> a.dividirPorEscalar(BigComplex.ZERO));
        assertThrows(DivisionPorCeroException.class, () -> a.cocienteHadamard(cero));
    }

    @Test void productoTrasposicionPermutacionYForma() {
        array.Matriz3D a = new array.Matriz3D(new double[][][]{
                {{1, 2, 3}, {4, 5, 6}},
                {{1, 0, 0}, {0, 1, 0}}
        });
        array.Matriz3D b = new array.Matriz3D(new double[][][]{
                {{1, 0}, {0, 1}, {1, 1}},
                {{2, 0}, {0, 2}, {0, 0}}
        });
        array.Matriz3D producto = a.productoPorCapas(b);
        assertEquals(new array.Matriz(new double[][]{{4, 5}, {10, 11}}), producto.capa(0));
        assertEquals(new array.Matriz(new double[][]{{2, 0}, {0, 2}}), producto.capa(1));
        assertEquals(producto, a.productoPorCapas(b.capa(0)).conCapa(1, producto.capa(1)));
        assertEquals(a.capa(0).traspuesta(), a.traspuestaPorCapas().capa(0));
        assertEquals(a.traspuestaPorCapas(), a.traspuesta());
        assertEquals("2x3x2", a.traspuestaPorCapas().dimensiones());

        array.Matriz3D permutada = a.permutarEjes(1, 0, 2);
        assertEquals("2x2x3", permutada.dimensiones());
        assertEquals(a.get(1, 0, 2), permutada.get(0, 1, 2));
        assertEquals(a, permutada.permutarEjes(1, 0, 2));
        assertThrows(IllegalArgumentException.class, () -> a.permutarEjes(0, 0, 2));

        array.Matriz3D reformada = a.reshape(3, 2, 2);
        assertArrayEquals(a.flatten(), reformada.flatten());
        assertEquals(a, reformada.reshape(2, 2, 3));
        assertTrue(a.serializar().startsWith("2x2x3:"));
        assertThrows(DimensionMatrizException.class, () -> a.reshape(2, 2, 2));
    }

    @Test void concatenacionFactoriasYApiBidimensionalPorCapa() {
        array.Matriz3D a = new array.Matriz3D(DATOS);
        assertEquals("4x2x3", a.concatenarCapas(a).dimensiones());
        assertEquals("2x4x3", a.concatenarFilas(a).dimensiones());
        assertEquals("2x2x6", a.concatenarColumnas(a).dimensiones());
        assertEquals(a, a.redimensionar(3, 3, 4, BigComplex.ZERO).submatriz(0, 2, 0, 2, 0, 3));

        array.Matriz3D identidades = array.Matriz3D.identidadPorCapas(3, 2, 2);
        for (int k = 0; k < identidades.capas(); k++) assertTrue(identidades.capa(k).esIdentidad());
        assertEquals(a.capa(0), array.Matriz3D.desdeMatriz(a.capa(0)).capa(0));
        assertEquals(a, array.Matriz3D.desdeCapas(a.capa(0), a.capa(1)));
        assertEquals(lista.Matriz3D.desdeCapas(a.capa(0), a.capa(1)), new lista.Matriz3D(a));
        assertEquals(array.Matriz3D.aleatoria(2, 2, 2, 7), array.Matriz3D.aleatoria(2, 2, 2, 7));

        array.Matriz3D cuadradas = new array.Matriz3D(new double[][][]{
                {{2, 0}, {0, 4}}, {{1, 0}, {0, 5}}
        });
        assertTrue(cuadradas.productoPorCapas(cuadradas.aplicarPorCapas(array.Matriz::inversa))
                .aproximadamenteIgual(array.Matriz3D.identidadPorCapas(2, 2, 2), TOL));
        assertEquals(cuadradas.suma(cuadradas), cuadradas.combinarPorCapas(cuadradas, array.Matriz::suma));
    }
}
