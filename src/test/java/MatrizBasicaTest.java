import ch.obermuhlner.math.big.BigComplex;
import comun.Complejos;
import comun.excepciones.DimensionMatrizException;
import comun.excepciones.DivisionPorCeroException;
import comun.excepciones.IndiceMatrizException;
import comun.excepciones.MatrizNoCuadradaException;
import comun.excepciones.MatrizSingularException;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.List;
import java.lang.reflect.Field;

import static org.junit.jupiter.api.Assertions.*;

class MatrizBasicaTest {
    private static final BigDecimal TOL = new BigDecimal("1E-20");

    @ParameterizedTest
    @ValueSource(ints = {1, 2, 4})
    void identidadesDeDistintosOrdenes(int orden) {
        array.Matriz a = array.Matriz.identidad(orden);
        lista.Matriz l = lista.Matriz.identidad(orden);
        assertTrue(a.esIdentidad());
        assertTrue(l.esIdentidad());
        assertEquals(a, l);
        assertEquals(BigComplex.ONE, a.determinante());
    }

    @Test void constructoresAccesoNoAliasingYConversion() {
        double[][] origen = {{1, 2}, {3, 4}};
        array.Matriz a = new array.Matriz(origen);
        origen[0][0] = 99;
        assertEquals(BigComplex.ONE, a.get(0, 0));
        BigComplex[][] copia = a.getMatriz();
        copia[0][0] = BigComplex.ZERO;
        assertEquals(BigComplex.ONE, a.get(0, 0));
        assertArrayEquals(new BigComplex[]{BigComplex.ONE, Complejos.de(2)}, a.obtenerFila(0));
        assertArrayEquals(new BigComplex[]{Complejos.de(2), Complejos.de(4)}, a.obtenerColumna(1));
        assertEquals(a, new array.Matriz(a));
        assertEquals(a, new array.Matriz(new Double[][]{{1.0, 2.0}, {3.0, 4.0}}));
        assertEquals(a, new array.Matriz(new BigDecimal[][]{{BigDecimal.ONE, BigDecimal.TWO},
                {BigDecimal.valueOf(3), BigDecimal.valueOf(4)}}));
        assertEquals(a, new array.Matriz(new String[][]{{"1", "2"}, {"3", "4"}}));
        assertTrue(new array.Matriz(2, 3, true).esIdentidadRectangular());
        assertTrue(new array.Matriz(2, 3).esNula());
        assertEquals(7, new array.Matriz(2, 2, Complejos.de(7)).get(0, 1).re.intValueExact());

        ArrayList<ArrayList<Double>> lista = new ArrayList<>();
        lista.add(new ArrayList<>(List.of(1.0, 2.0)));
        lista.add(new ArrayList<>(List.of(3.0, 4.0)));
        lista.Matriz l = new lista.Matriz(lista);
        lista.getFirst().set(0, 50.0);
        assertEquals(a, l);
        assertEquals(a, l.comoMatrizArray());
        assertEquals(l, a.comoMatrizLista());
        assertNotSame(l.comoListaProfunda(), l.comoListaProfunda());
    }

    @Test void cadaImplementacionConservaSuAlmacenamientoReal() throws Exception {
        Field campoArray = array.Matriz.class.getDeclaredField("matriz");
        Field campoLista = lista.Matriz.class.getDeclaredField("matriz");
        assertEquals(BigComplex[][].class, campoArray.getType());
        assertEquals(ArrayList.class, campoLista.getType());
        lista.Matriz a = new lista.Matriz(new double[][]{{1, 2}, {3, 4}});
        lista.Matriz resultado = a.suma(a).producto(a).traspuesta().inversa();
        assertInstanceOf(lista.Matriz.class, resultado);
        assertFalse(java.nio.file.Files.readString(java.nio.file.Path.of(
                "src/main/java/lista/Matriz.java")).contains("array.Matriz"));
    }

    @Test void accesoFuncionalYDimensiones() {
        array.Matriz a = new array.Matriz(new double[][]{{1, 2}, {3, 4}});
        assertEquals(9, a.conElemento(0, 0, Complejos.de(9)).get(0, 0).re.intValue());
        assertEquals(8, a.conFila(0, new BigComplex[]{Complejos.de(8), Complejos.de(7)}).get(0, 0).re.intValue());
        assertEquals(6, a.conColumna(1, new BigComplex[]{Complejos.de(6), Complejos.de(5)}).get(0, 1).re.intValue());
        assertEquals(1, a.submatriz(0, 1, 0, 1).getFilas());
        assertEquals(1, a.eliminarFila(0).getFilas());
        assertEquals(1, a.eliminarColumna(0).getColumnas());
        assertEquals(4, a.concatenarHorizontal(a).getColumnas());
        assertEquals(4, a.concatenarVertical(a).getFilas());
        assertEquals(3, a.insertarFila(1, a.obtenerFila(0)).getFilas());
        assertEquals(3, a.insertarColumna(1, a.obtenerColumna(0)).getColumnas());
        assertEquals(3, a.redimensionar(3, 3, BigComplex.ZERO).getFilas());
        assertTrue(a.rellenar(BigComplex.ZERO).esNula());
        assertEquals(a.productoEscalar(2), a.map(z -> z.multiply(2)));
        assertEquals(a.suma(a), a.zip(a, BigComplex::add));
        assertThrows(IndiceMatrizException.class, () -> a.get(-1, 0));
        assertThrows(IndiceMatrizException.class, () -> a.get(2, 0));
        assertThrows(DimensionMatrizException.class, () -> new array.Matriz(new double[][]{{1}, {2, 3}}));
        assertThrows(NullPointerException.class, () -> new array.Matriz(new Double[][]{{1.0, null}}));
    }

    @Test void aritmeticaYPropiedades() {
        array.Matriz a = new array.Matriz(new double[][]{{1, 2}, {3, 4}});
        array.Matriz b = new array.Matriz(new double[][]{{0, 1}, {1, 0}});
        array.Matriz cero = array.Matriz.cero(2, 2);
        assertEquals(a, a.suma(cero));
        assertTrue(a.resta(a).esNula());
        assertTrue(a.productoEscalar(BigComplex.ZERO).esNula());
        assertEquals(a, a.productoEscalar(BigComplex.ONE));
        assertEquals(a, a.dividirPorEscalar(BigDecimal.ONE));
        assertEquals(a.producto(b.suma(a)), a.producto(b).suma(a.producto(a)));
        assertEquals(a.producto(b).producto(a), a.producto(b.producto(a)));
        assertNotEquals(a.producto(b), b.producto(a));
        assertEquals(a.producto(b).traspuesta(), b.traspuesta().producto(a.traspuesta()));
        assertEquals(a.producto(b).adjuntaHermitica(), b.adjuntaHermitica().producto(a.adjuntaHermitica()));
        assertEquals(a.producto(a), a.potencia(2));
        assertTrue(a.potencia(0).esIdentidad());
        assertTrue(a.potencia(-1).producto(a).aproximadamenteIgual(array.Matriz.identidad(2), TOL));
        assertEquals(a.productoHadamard(b), new array.Matriz(new double[][]{{0, 2}, {3, 0}}));
        assertEquals(4, a.productoKronecker(b).getFilas());
        assertEquals(4, a.productoKhatriRao(b).getFilas());
        assertThrows(DivisionPorCeroException.class, () -> a.dividirPorEscalar(BigComplex.ZERO));
        assertThrows(DivisionPorCeroException.class, () -> a.cocienteHadamard(cero));
        assertThrows(DimensionMatrizException.class, () -> a.suma(new array.Matriz(1, 2)));
        assertThrows(DimensionMatrizException.class, () -> a.producto(new array.Matriz(3, 1)));
    }

    @Test void transposicionPredicadosYEspeciales() {
        array.Matriz h = new array.Matriz(new String[][]{{"2", "1+i"}, {"1-i", "3"}});
        assertTrue(h.esHermitica());
        assertEquals(h, h.adjuntaHermitica());
        assertTrue(h.parteAntihermitica().esNula());
        assertTrue(h.parteHermitica().aproximadamenteIgual(h, TOL));
        assertTrue(array.Matriz.identidad(3).esDiagonal());
        assertTrue(array.Matriz.identidad(3).esEscalar());
        assertTrue(array.Matriz.hadamard(4).esMatrizDeSignos());
        assertTrue(array.Matriz.matrizDePermutacion(1, 0).esMatrizDePermutacion());
        assertTrue(array.Matriz.toeplitz(new BigComplex[]{Complejos.de(1), Complejos.de(2)},
                new BigComplex[]{Complejos.de(1), Complejos.de(3)}).esToeplitz());
        assertTrue(array.Matriz.hankel(new BigComplex[]{Complejos.de(1), Complejos.de(2)},
                new BigComplex[]{Complejos.de(2), Complejos.de(3)}).esHankel());
        assertTrue(array.Matriz.circulante(Complejos.de(1), Complejos.de(2)).esCirculante());
        assertTrue(array.Matriz.matrizDeRotacion2D(BigDecimal.ZERO).esRotacion());
        assertTrue(array.Matriz.matrizDeReflexion(new BigComplex[]{BigComplex.ONE, BigComplex.ZERO}).esReflexion());
        assertTrue(array.Matriz.matrizDeProyeccion(new BigComplex[]{BigComplex.ONE, BigComplex.ZERO}).esProyectorOrtogonal());
    }

    @Test void gaussDeterminanteInversaRangoYSistemas() {
        array.Matriz a = new array.Matriz(new double[][]{{4, 7}, {2, 6}});
        assertEquals(Complejos.de(10), a.determinante());
        assertEquals(a.determinante(), a.determinanteLaplace());
        assertEquals(a.determinante(), a.determinanteLU());
        assertEquals(2, a.rango());
        assertEquals(0, a.nulidad());
        assertTrue(a.verificarTeoremaRangoNulidad());
        assertTrue(a.rref().esIdentidad());
        assertTrue(a.verificarInversa(a.inversa(), TOL));
        assertTrue(a.verificarInversa(a.inversaPorAdjunta(), TOL));
        array.Matriz b = array.Matriz.columna(Complejos.de(1), Complejos.de(0));
        array.Matriz x = a.resolver(b);
        assertTrue(a.verificarSolucionSistema(x, b, TOL));
        assertEquals(comun.resultados.TipoSolucion.UNICA, a.discutirSistema(b).tipo());
        assertEquals(a.rango(), a.columnasPivote().size());
        assertThrows(MatrizNoCuadradaException.class, () -> new array.Matriz(2, 3).determinante());
        assertThrows(MatrizSingularException.class, () -> new array.Matriz(new double[][]{{1, 2}, {2, 4}}).inversa());
    }
}
