import ch.obermuhlner.math.big.BigComplex;
import comun.Complejos;
import comun.ModoConvolucion;
import comun.excepciones.DimensionMatrizException;
import comun.resultados.ResultadoHOSVD;
import comun.resultados.ResultadoTSVD;
import org.junit.jupiter.api.Test;

import java.math.BigDecimal;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

class Matriz3DOperacionesTest {
    private static final double[][][] DATOS = {
            {{1, 2, 3}, {4, 5, 6}},
            {{7, 8, 9}, {10, 11, 12}}
    };

    @Test void fibrasPlanosEIndicesSobreTodosLosEjes() {
        array.Matriz3D a = new array.Matriz3D(DATOS);
        assertArrayEquals(z(2, 8), a.fibra(0, 0, 1));
        assertArrayEquals(z(9, 12), a.fibra(1, 1, 2));
        assertArrayEquals(z(4, 5, 6), a.fibra(2, 0, 1));
        assertEquals(a.capa(1), a.plano(0, 1));
        assertEquals(a.cortePorFila(1), a.plano(1, 1));
        assertEquals(a.cortePorColumna(2), a.plano(2, 2));

        array.Matriz3D editada = a.conFibra(2, 0, 0, z(30, 20, 10));
        assertArrayEquals(z(30, 20, 10), editada.fibra(2, 0, 0));
        assertArrayEquals(z(12, 9), a.extraerIndices(1, 1, 0).fibra(1, 1, 2));
        assertEquals(a, a.insertarPlano(0, 1, a.capa(0)).eliminarIndice(0, 1));
        assertEquals(a, a.insertarPlano(1, 1, a.cortePorFila(0)).eliminarIndice(1, 1));
        assertEquals(a, a.insertarPlano(2, 1, a.cortePorColumna(0)).eliminarIndice(2, 1));

        List<array.Matriz3D> partes = a.particionar(2, 1, 2);
        assertEquals(3, partes.size());
        assertEquals(a, partes.get(0).concatenarColumnas(partes.get(1)).concatenarColumnas(partes.get(2)));
    }

    @Test void broadcastingPaddingRepeticionYReordenacion() {
        array.Matriz3D fila = new array.Matriz3D(new double[][][]{{{1, 2}}});
        array.Matriz3D expandida = fila.broadcast(2, 3, 2);
        assertEquals("2x3x2", expandida.dimensiones());
        assertArrayEquals(z(1, 2), expandida.fibra(2, 1, 2));
        assertEquals(expandida, fila.sumaBroadcast(array.Matriz3D.cero(2, 3, 2)));
        assertEquals(fila.repetir(2, 3, 1), expandida);

        array.Matriz3D borde = fila.rellenarBordes(1, 1, 1, 0, 2, 1, Complejos.de(-1));
        assertEquals("3x2x5", borde.dimensiones());
        assertEquals(BigComplex.ONE, borde.get(1, 1, 2));
        assertEquals(Complejos.de(-1), borde.get(0, 0, 0));

        array.Matriz3D a = new array.Matriz3D(DATOS);
        assertArrayEquals(z(9, 8, 7), a.voltear(2).fibra(2, 1, 0));
        assertArrayEquals(z(3, 1, 2), new array.Matriz3D(new double[][][]{{{1, 2, 3}}})
                .desplazarCiclicamente(1, 2).fibra(2, 0, 0));
        assertEquals(a, a.moverEje(2, 0).moverEje(0, 2));
        assertEquals(new array.Matriz(new double[][]{{3, 6}, {2, 5}, {1, 4}}),
                a.extraerCapas(0).rotar90(1, 1, 2).capa(0));
        assertEquals(a, a.rotar90(4, 0, 1));
    }

    @Test void reduccionesEstadisticasYOperacionesElemento() {
        array.Matriz3D a = new array.Matriz3D(new double[][][]{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}});
        assertEquals(Complejos.de(40320), a.productoElementos());
        assertEquals(new array.Matriz(new double[][]{{6, 8}, {10, 12}}), a.sumaEje(0));
        assertEquals(new array.Matriz(new double[][]{{4, 6}, {12, 14}}), a.sumaEje(1));
        assertEquals(new array.Matriz(new double[][]{{3, 7}, {11, 15}}), a.sumaEje(2));
        assertEquals(a.sumaPorCapas(), a.sumaEje(0));
        assertEquals(new array.Matriz(new double[][]{{3, 4}, {5, 6}}),
                a.promedioPonderadoEje(0, new BigDecimal[]{BigDecimal.ONE, BigDecimal.ONE}));
        array.Matriz normas = a.normaFrobeniusEje(0);
        assertTrue(a.productoInternoEje(a, 0).aproximadamenteIgual(normas.productoHadamard(normas),
                new BigDecimal("1E-20")));
        assertEquals(BigDecimal.valueOf(5.25), a.varianza());
        assertEquals(0, BigDecimal.valueOf(8).compareTo(a.elementoMayorModulo().re));
        assertEquals(BigComplex.ONE, a.elementoMenorModulo());
        assertArrayEquals(new int[]{1, 1, 1}, a.indiceMayorModulo());
        assertEquals(8, a.contarNoNulos());
        assertEquals(a.productoHadamard(a), a.potenciaElementos(2));
        assertEquals(a, a.parteReal());
        assertTrue(a.parteImaginaria().esNula());
        assertEquals(a, a.valorAbsoluto());
    }

    @Test void algebraTensorialModoContraccionesYProductos() {
        array.Matriz3D a = new array.Matriz3D(new double[][][]{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}});
        assertEquals("2x4", a.desplegar(0).filas() + "x" + a.desplegar(0).columnas());
        assertArrayEquals(z(1, 2, 3, 4), a.desplegar(0).obtenerFila(0));
        array.Matriz suma = new array.Matriz(new double[][]{{1, 1}});
        array.Matriz3D modo0 = a.productoModo(suma, 0);
        assertEquals("1x2x2", modo0.dimensiones());
        assertEquals(new array.Matriz(new double[][]{{6, 8}, {10, 12}}), modo0.capa(0));
        assertEquals(a.sumaEje(0), a.contraerModo(z(1, 1), 0));
        assertEquals(Complejos.de(36), a.formaMultilineal(z(1, 1), z(1, 1), z(1, 1)));
        assertArrayEquals(new int[]{2, 2, 2}, a.rangoMultilineal());
        assertArrayEquals(z(5, 13), a.traza(1, 2));

        array.Matriz3D diagonal = array.Matriz3D.cero(2, 2, 2)
                .conElemento(0, 0, 0, BigComplex.ONE).conElemento(1, 1, 1, Complejos.de(2));
        assertTrue(diagonal.esDiagonal());
        assertTrue(diagonal.esSimetrica());
        assertArrayEquals(z(1, 2), diagonal.diagonalPrincipal());
        assertEquals(diagonal, diagonal.simetrizar());
        assertTrue(diagonal.antisimetrizar().esNula());

        array.Matriz3D unidad = new array.Matriz3D(new double[][][]{{{2}}});
        assertEquals(a.productoEscalar(2), a.productoKronecker(unidad));
        assertEquals(Complejos.de(30), a.contraerDoble(a, 1, 2, 1, 2).get(0, 0));
        assertEquals(Complejos.de(16), a.productoExterior(unidad)[1][1][1][0][0][0]);

        array.Matriz3D x = new array.Matriz3D(new double[][][]{{{1, 0, 0}}});
        array.Matriz3D y = new array.Matriz3D(new double[][][]{{{0, 1, 0}}});
        assertArrayEquals(z(0, 0, 1), x.productoVectorial(y, 2).fibra(2, 0, 0));
    }

    @Test void convolucionYCorrelacionTridimensional() {
        array.Matriz3D senal = new array.Matriz3D(new double[][][]{{{1, 2, 3}}});
        array.Matriz3D nucleo = new array.Matriz3D(new double[][][]{{{1, 1}}});
        assertArrayEquals(z(1, 3, 5, 3), senal.convolucion(nucleo, ModoConvolucion.COMPLETA).fibra(2, 0, 0));
        assertArrayEquals(z(1, 3, 5), senal.convolucion(nucleo, ModoConvolucion.MISMA).fibra(2, 0, 0));
        assertArrayEquals(z(3, 5), senal.convolucion(nucleo, ModoConvolucion.VALIDA).fibra(2, 0, 0));
        assertEquals(senal.convolucion(nucleo, ModoConvolucion.COMPLETA),
                senal.correlacion(nucleo, ModoConvolucion.COMPLETA));
        assertThrows(DimensionMatrizException.class, () -> nucleo.convolucion(senal, ModoConvolucion.VALIDA));
        assertInstanceOf(lista.Matriz3D.class,
                new lista.Matriz3D(senal).convolucion(new lista.Matriz3D(nucleo), ModoConvolucion.MISMA));
    }

    @Test void factoriasSerializacionYHosvd() {
        array.Matriz3D secuencia = array.Matriz3D.secuencial(2, 2, 2, BigComplex.ONE, BigComplex.ONE);
        assertArrayEquals(z(1, 2, 3, 4, 5, 6, 7, 8), secuencia.flatten());
        assertEquals(secuencia, array.Matriz3D.desdeVector(secuencia.flatten(), 2, 2, 2));
        assertEquals(secuencia, array.Matriz3D.deserializar(secuencia.serializar()));
        assertEquals(new lista.Matriz3D(secuencia), lista.Matriz3D.deserializar(secuencia.serializar()));
        assertTrue(array.Matriz3D.unos(2, 2, 2).resta(array.Matriz3D.rellena(2, 2, 2, BigComplex.ONE)).esNula());
        assertTrue(array.Matriz3D.identidad(3).esDiagonal());
        assertArrayEquals(z(2, 3), array.Matriz3D.diagonal(z(2, 3)).diagonalPrincipal());
        assertEquals(array.Matriz3D.aleatoriaCompleja(2, 2, 2, 4),
                array.Matriz3D.aleatoriaCompleja(2, 2, 2, 4));
        assertEquals(0, array.Matriz3D.aleatoriaDispersa(2, 2, 2, BigDecimal.ZERO, 4).contarNoNulos());

        array.Matriz3D diagonal = array.Matriz3D.diagonal(z(2, 3));
        ResultadoHOSVD<array.Matriz3D, array.Matriz> hosvd = diagonal.descomposicionHOSVD();
        assertTrue(diagonal.verificarHOSVD(hosvd, new BigDecimal("1E-15")));
        assertEquals("1x1x1", diagonal.descomposicionHOSVDTruncada(1, 1, 1).nucleo().dimensiones());
    }

    @Test void algebraTProduct() {
        array.Matriz3D a = new array.Matriz3D(new double[][][]{{{1}}, {{2}}});
        array.Matriz3D b = new array.Matriz3D(new double[][][]{{{3}}, {{4}}});
        assertArrayEquals(z(11, 10), a.productoT(b).fibra(0, 0, 0));
        array.Matriz3D identidad = array.Matriz3D.identidadT(2, 1);
        assertTrue(identidad.esIdentidadT());
        assertEquals(a, a.productoT(identidad));
        assertTrue(a.transformadaFourierCapas().transformadaInversaFourierCapas()
                .aproximadamenteIgual(a, new BigDecimal("1E-20")));
        array.Matriz3D escalar = new array.Matriz3D(new double[][][]{{{2}}, {{0}}});
        assertTrue(escalar.productoT(escalar.inversaT()).aproximadamenteIgual(identidad, new BigDecimal("1E-20")));
        assertTrue(escalar.potenciaT(0).esIdentidadT());
        assertTrue(escalar.potenciaT(2).aproximadamenteIgual(escalar.productoT(escalar), new BigDecimal("1E-20")));
        assertEquals(a.adjuntaT().adjuntaT(), a);
        ResultadoTSVD<array.Matriz3D> tsvd = a.descomposicionTSVD();
        assertTrue(a.verificarTSVD(tsvd, new BigDecimal("1E-20")));
    }

    private static BigComplex[] z(double... valores) {
        BigComplex[] salida = new BigComplex[valores.length];
        for (int i = 0; i < valores.length; i++) salida[i] = Complejos.de(valores[i]);
        return salida;
    }
}
