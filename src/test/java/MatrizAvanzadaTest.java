import ch.obermuhlner.math.big.BigComplex;
import comun.Complejos;
import comun.resultados.ResultadoLU;
import comun.resultados.ResultadoQR;
import org.junit.jupiter.api.Test;

import java.math.BigDecimal;

import static org.junit.jupiter.api.Assertions.*;

class MatrizAvanzadaTest {
    private static final BigDecimal TOL = new BigDecimal("1E-18");

    @Test void factorizacionesYNormas() {
        array.Matriz a = new array.Matriz(new double[][]{{4, 1}, {1, 3}});
        ResultadoLU<array.Matriz> lu = a.descomposicionLUP();
        assertTrue(a.verificarLU(lu, TOL));
        ResultadoQR<array.Matriz> qr = a.descomposicionQRGramSchmidtModificado();
        assertTrue(a.verificarQR(qr, TOL));
        ResultadoLU<array.Matriz> ch = a.descomposicionCholesky();
        assertTrue(a.verificarCholesky(ch, TOL));
        assertTrue(a.esDefinidaPositiva());
        assertTrue(a.normaFrobenius().signum() > 0);
        assertTrue(a.normaUno().signum() > 0);
        assertTrue(a.normaInfinito().signum() > 0);
        assertEquals(a.traza(), a.traspuesta().traza());
        assertEquals(a.trazaDePotencia(2), a.producto(a).traza());
    }

    @Test void diagonalEspectroFuncionesYCalculo() {
        array.Matriz d = array.Matriz.diagonal(Complejos.de(2), Complejos.de(3));
        assertTrue(d.esDiagonalizable());
        var diag = d.diagonalizacion();
        assertTrue(d.verificarDiagonalizacion(diag, TOL));
        assertTrue(d.verificarCayleyHamilton(TOL));
        assertEquals(2, d.autovalores().size());
        assertTrue(d.exponencial().aproximadamenteIgual(d.exponencialPorDiagonalizacion(), new BigDecimal("1E-15")));
        assertTrue(d.raizCuadradaPrincipal().potencia(2).aproximadamenteIgual(d, TOL));
        assertTrue(d.seno().producto(d.seno()).suma(d.coseno().producto(d.coseno()))
                .aproximadamenteIgual(array.Matriz.identidad(2), new BigDecimal("1E-15")));
        assertEquals(d, d.derivadaNumerica(t -> d.productoEscalar(t), BigDecimal.ONE,
                new BigDecimal("1E-10")).map(z -> z.round(d.mathContext())));
    }

    @Test void iterativosEstadisticaGrafosMarkovVectorizacion() {
        array.Matriz a = new array.Matriz(new double[][]{{4, 1}, {1, 3}});
        array.Matriz b = array.Matriz.columna(Complejos.de(1), Complejos.de(2));
        assertTrue(a.resolverJacobi(b, TOL, 200).convergio());
        assertTrue(a.resolverGaussSeidel(b, TOL, 200).convergio());
        assertTrue(a.gradienteConjugado(b, TOL, 20).convergio());
        array.Matriz datos = new array.Matriz(new double[][]{{1, 2}, {2, 4}, {3, 8}});
        assertTrue(datos.centrarColumnas().obtenerColumna(0)[0].re.signum() < 0);
        assertTrue(datos.matrizCovarianza().esHermitica());
        array.Matriz ady = new array.Matriz(new double[][]{{0, 1}, {1, 0}});
        assertTrue(ady.esGrafoConexo());
        assertEquals(1, ady.componentesConexas().size());
        assertTrue(ady.laplaciana().esSemidefinidaPositiva());
        array.Matriz p = new array.Matriz(new double[][]{{0.5, 0.5}, {0.25, 0.75}});
        assertTrue(p.esMatrizTransicion());
        assertTrue(p.esCadenaRegular());
        assertEquals(2, p.distribucionEstacionaria().getColumnas());
        assertEquals(a, a.vectorizarPorFilas().desvectorizar(2, 2, false));
        assertEquals(a, a.reshape(2, 2));
    }
}
