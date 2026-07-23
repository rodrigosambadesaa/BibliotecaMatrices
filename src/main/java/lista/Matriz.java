package lista;

import ch.obermuhlner.math.big.BigComplex;
import comun.AbstractMatriz;
import comun.Complejos;
import comun.MatrizCompleja;
import comun.excepciones.DimensionMatrizException;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;
import java.util.Random;

/**
 * Matriz inmutable almacenada realmente en
 * {@code ArrayList<ArrayList<BigComplex>>}. Todas las estructuras se copian,
 * nunca se exponen, y comparte únicamente algoritmos de acceso con la versión
 * array. La política predeterminada es DECIMAL128 con tolerancia 1E-28.
 */
public final class Matriz extends AbstractMatriz<Matriz> {
    private final ArrayList<ArrayList<BigComplex>> matriz;

    public Matriz(int filas, int columnas, boolean diagonal) {
        this(datosDimensiones(filas, columnas, BigComplex.ZERO, diagonal));
    }
    public Matriz(double[][] array) { this(convertir(array)); }
    public Matriz(Double[][] array) { this(convertir(array)); }
    public Matriz(BigDecimal[][] array) { this(convertir(array)); }
    public Matriz(BigComplex[][] array) { this(array, Complejos.CONTEXTO_DEFECTO, Complejos.TOLERANCIA_DEFECTO); }
    public Matriz(String[][] array) { this(convertir(array)); }
    public Matriz(Matriz m) { this((MatrizCompleja) Objects.requireNonNull(m)); }
    public Matriz(int filas, int columnas) { this(filas, columnas, false); }
    public Matriz(int filas, int columnas, BigComplex valorInicial) {
        this(datosDimensiones(filas, columnas, valorInicial, false));
    }
    public Matriz(int filas, int columnas, MathContext mathContext) {
        this(datosDimensiones(filas, columnas, BigComplex.ZERO, false), mathContext,
                BigDecimal.ONE.scaleByPowerOfTen(-Math.max(1, mathContext.getPrecision() - 6)));
    }
    public Matriz(ArrayList<ArrayList<Double>> v) { this(convertirDouble(v)); }
    public Matriz(List<? extends List<? extends Number>> v) { this(convertirNumber(v)); }
    public Matriz(MatrizCompleja m) { this(m, m.mathContext(), m.tolerancia()); }
    public Matriz(MatrizCompleja m, MathContext contexto, BigDecimal tolerancia) {
        this(copiar(m), contexto, tolerancia);
    }
    public Matriz(BigComplex[][] datos, MathContext contexto, BigDecimal tolerancia) {
        this(validar(datos, contexto), contexto, tolerancia, true);
    }
    private Matriz(BigComplex[][] datos, MathContext contexto, BigDecimal tolerancia, boolean validados) {
        super(datos.length, datos[0].length, contexto, tolerancia);
        matriz = new ArrayList<>(datos.length);
        for (BigComplex[] fila : datos) {
            ArrayList<BigComplex> f = new ArrayList<>(fila.length);
            for (BigComplex z : fila) f.add(validados ? z : z.round(contexto));
            matriz.add(f);
        }
    }

    public static Matriz desdeBigComplex(ArrayList<ArrayList<BigComplex>> v) {
        Objects.requireNonNull(v);
        if (v.isEmpty() || v.getFirst() == null || v.getFirst().isEmpty()) throw new DimensionMatrizException("Lista vacía");
        BigComplex[][] d = new BigComplex[v.size()][v.getFirst().size()];
        for (int i = 0; i < v.size(); i++) {
            if (v.get(i) == null || v.get(i).size() != d[0].length) throw new DimensionMatrizException("Lista irregular");
            for (int j = 0; j < d[i].length; j++) d[i][j] = Objects.requireNonNull(v.get(i).get(j));
        }
        return new Matriz(d);
    }

    @Override public BigComplex get(int fila, int columna) {
        validarIndice(fila, columna); return matriz.get(fila).get(columna);
    }
    @Override protected Matriz crear(BigComplex[][] datos, MathContext contexto, BigDecimal tol) {
        return new Matriz(datos, contexto, tol);
    }

    public static Matriz cero(int f, int c) { return new Matriz(f, c); }
    public static Matriz identidad(int n) { return new Matriz(n, n, true); }
    public static Matriz identidad(int f, int c) { return new Matriz(f, c, true); }
    public static Matriz diagonal(BigComplex... e) {
        if (e == null || e.length == 0) throw new DimensionMatrizException("Diagonal vacía");
        BigComplex[][] d = datosDimensiones(e.length, e.length, BigComplex.ZERO, false);
        for (int i = 0; i < e.length; i++) d[i][i] = Objects.requireNonNull(e[i]);
        return new Matriz(d);
    }
    public static Matriz escalar(int n, BigComplex v) {
        BigComplex[] d = new BigComplex[n]; Arrays.fill(d, Objects.requireNonNull(v)); return diagonal(d);
    }
    public static Matriz fila(BigComplex... v) { return new Matriz(new BigComplex[][]{v}); }
    public static Matriz columna(BigComplex... v) {
        if (v == null || v.length == 0) throw new DimensionMatrizException("Columna vacía");
        BigComplex[][] d = new BigComplex[v.length][1];
        for (int i = 0; i < v.length; i++) d[i][0] = v[i];
        return new Matriz(d);
    }
    public static Matriz desdeFilas(BigComplex[]... v) { return new Matriz(v); }
    public static Matriz desdeColumnas(BigComplex[]... v) {
        if (v == null || v.length == 0 || v[0] == null) throw new DimensionMatrizException("Columnas vacías");
        BigComplex[][] d = new BigComplex[v[0].length][v.length];
        for (int j = 0; j < v.length; j++) {
            if (v[j] == null || v[j].length != d.length) throw new DimensionMatrizException("Columnas irregulares");
            for (int i = 0; i < d.length; i++) d[i][j] = v[j][i];
        }
        return new Matriz(d);
    }
    public static Matriz aleatoria(int f, int c, long s) {
        Random r = new Random(s); BigComplex[][] d = new BigComplex[f][c];
        for (int i = 0; i < f; i++) for (int j = 0; j < c; j++)
            d[i][j] = BigComplex.valueOf(BigDecimal.valueOf(r.nextDouble() * 2 - 1));
        return new Matriz(d);
    }
    public static Matriz aleatoriaHermitica(int n, long s) {
        Matriz a = aleatoria(n, n, s); return a.suma(a.adjuntaHermitica()).dividirPorEscalar(BigDecimal.TWO);
    }
    public static Matriz aleatoriaUnitaria(int n, long s) {
        return aleatoria(n, n, s).descomposicionQRGramSchmidtModificado().q();
    }
    public static Matriz aleatoriaDefinidaPositiva(int n, long s) {
        Matriz a = aleatoria(n, n, s); return a.adjuntaHermitica().producto(a).suma(identidad(n));
    }
    public static Matriz aleatoriaDispersa(int f, int c, BigDecimal densidad, long s) {
        if (densidad.signum() < 0 || densidad.compareTo(BigDecimal.ONE) > 0) throw new IllegalArgumentException("Densidad inválida");
        Random r = new Random(s); BigComplex[][] d = datosDimensiones(f, c, BigComplex.ZERO, false);
        for (int i = 0; i < f; i++) for (int j = 0; j < c; j++)
            if (r.nextDouble() < densidad.doubleValue()) d[i][j] = Complejos.de(r.nextDouble());
        return new Matriz(d);
    }
    public static Matriz matrizDePermutacion(int... p) {
        if (p == null || Arrays.stream(p).distinct().count() != p.length
                || Arrays.stream(p).anyMatch(x -> x < 0 || x >= p.length)) throw new DimensionMatrizException("Permutación inválida");
        BigComplex[][] d = datosDimensiones(p.length, p.length, BigComplex.ZERO, false);
        for (int i = 0; i < p.length; i++) d[i][p[i]] = BigComplex.ONE;
        return new Matriz(d);
    }
    public static Matriz matrizDeRotacion2D(BigDecimal a) {
        MathContext mc = Complejos.CONTEXTO_DEFECTO;
        BigDecimal c = new BigDecimal(Math.cos(a.doubleValue()), mc), s = new BigDecimal(Math.sin(a.doubleValue()), mc);
        return new Matriz(new BigDecimal[][]{{c, s.negate()}, {s, c}});
    }
    public static Matriz matrizDeRotacion3D(BigDecimal a, char e) {
        Matriz r = matrizDeRotacion2D(a); BigComplex c = r.get(0, 0), s = r.get(1, 0);
        return switch (Character.toLowerCase(e)) {
            case 'x' -> new Matriz(new BigComplex[][]{{BigComplex.ONE, BigComplex.ZERO, BigComplex.ZERO},
                    {BigComplex.ZERO, c, s.negate()}, {BigComplex.ZERO, s, c}});
            case 'y' -> new Matriz(new BigComplex[][]{{c, BigComplex.ZERO, s}, {BigComplex.ZERO, BigComplex.ONE, BigComplex.ZERO},
                    {s.negate(), BigComplex.ZERO, c}});
            case 'z' -> new Matriz(new BigComplex[][]{{c, s.negate(), BigComplex.ZERO}, {s, c, BigComplex.ZERO},
                    {BigComplex.ZERO, BigComplex.ZERO, BigComplex.ONE}});
            default -> throw new IllegalArgumentException("Eje inválido");
        };
    }
    public static Matriz matrizDeReflexion(BigComplex[] n) {
        Matriz v = columna(n).normalizarColumnas();
        return identidad(n.length).resta(v.producto(v.adjuntaHermitica()).productoEscalar(BigDecimal.TWO));
    }
    public static Matriz matrizDeProyeccion(BigComplex[] d) {
        Matriz v = columna(d).normalizarColumnas(); return v.producto(v.adjuntaHermitica());
    }
    public static Matriz householder(BigComplex[] v) { return matrizDeReflexion(v); }
    public static Matriz givens(int n, int i, int j, BigDecimal a) {
        Matriz g = identidad(n), r = matrizDeRotacion2D(a);
        return g.conElemento(i, i, r.get(0, 0)).conElemento(i, j, r.get(0, 1))
                .conElemento(j, i, r.get(1, 0)).conElemento(j, j, r.get(1, 1));
    }
    public static Matriz toeplitz(BigComplex[] c, BigComplex[] f) {
        if (c == null || f == null || c.length == 0 || f.length == 0 || !Objects.equals(c[0], f[0]))
            throw new DimensionMatrizException("Datos Toeplitz inválidos");
        BigComplex[][] d = new BigComplex[c.length][f.length];
        for (int i = 0; i < d.length; i++) for (int j = 0; j < d[i].length; j++) d[i][j] = j >= i ? f[j - i] : c[i - j];
        return new Matriz(d);
    }
    public static Matriz hankel(BigComplex[] f, BigComplex[] c) {
        if (f == null || c == null || f.length == 0 || c.length == 0 || !Objects.equals(f[f.length - 1], c[0]))
            throw new DimensionMatrizException("Datos Hankel inválidos");
        BigComplex[] sec = new BigComplex[f.length + c.length - 1];
        System.arraycopy(f, 0, sec, 0, f.length); System.arraycopy(c, 1, sec, f.length, c.length - 1);
        BigComplex[][] d = new BigComplex[c.length][f.length];
        for (int i = 0; i < d.length; i++) for (int j = 0; j < d[i].length; j++) d[i][j] = sec[i + j];
        return new Matriz(d);
    }
    public static Matriz circulante(BigComplex... f) {
        int n = f.length; BigComplex[][] d = new BigComplex[n][n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) d[i][j] = f[(j - i + n) % n];
        return new Matriz(d);
    }
    public static Matriz vandermonde(BigComplex... nodos) {
        BigComplex[][] d = new BigComplex[nodos.length][nodos.length];
        for (int i = 0; i < nodos.length; i++) {
            d[i][0] = BigComplex.ONE;
            for (int j = 1; j < nodos.length; j++) d[i][j] = d[i][j - 1].multiply(nodos[i], Complejos.CONTEXTO_DEFECTO);
        }
        return new Matriz(d);
    }
    public static Matriz hilbert(int n) {
        BigDecimal[][] d = new BigDecimal[n][n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
            d[i][j] = BigDecimal.ONE.divide(BigDecimal.valueOf(i + j + 1L), Complejos.CONTEXTO_DEFECTO);
        return new Matriz(d);
    }
    public static Matriz fourierDiscreta(int n) {
        BigComplex[][] d = new BigComplex[n][n]; MathContext mc = Complejos.CONTEXTO_DEFECTO;
        BigDecimal escala = BigDecimal.ONE.divide(BigDecimal.valueOf(n).sqrt(mc), mc);
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
            d[i][j] = BigComplex.valueOfPolar(escala, new BigDecimal(-2 * Math.PI * i * j / n, mc), mc);
        return new Matriz(d);
    }
    public static Matriz hadamard(int n) {
        if (n <= 0 || (n & (n - 1)) != 0) throw new IllegalArgumentException("Sylvester exige potencia de dos");
        BigComplex[][] d = {{BigComplex.ONE}};
        while (d.length < n) {
            int m = d.length; BigComplex[][] x = new BigComplex[2 * m][2 * m];
            for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) {
                x[i][j] = d[i][j]; x[i][j + m] = d[i][j]; x[i + m][j] = d[i][j]; x[i + m][j + m] = d[i][j].negate();
            }
            d = x;
        }
        return new Matriz(d);
    }
    public static Matriz triangularSuperior(BigComplex[][] datos) {
        Matriz m = new Matriz(datos);
        if (!m.esTriangularSuperior()) throw new DimensionMatrizException("Los datos no son triangulares superiores");
        return m;
    }
    public static Matriz triangularInferior(BigComplex[][] datos) {
        Matriz m = new Matriz(datos);
        if (!m.esTriangularInferior()) throw new DimensionMatrizException("Los datos no son triangulares inferiores");
        return m;
    }
    public static Matriz bidiagonalSuperior(BigComplex[] diagonal, BigComplex[] superior) {
        if (diagonal == null || superior == null || superior.length != diagonal.length - 1)
            throw new DimensionMatrizException("Bandas incompatibles");
        BigComplex[][] d = datosDimensiones(diagonal.length, diagonal.length, BigComplex.ZERO, false);
        for (int i = 0; i < diagonal.length; i++) {
            d[i][i] = Objects.requireNonNull(diagonal[i]);
            if (i + 1 < diagonal.length) d[i][i + 1] = Objects.requireNonNull(superior[i]);
        }
        return new Matriz(d);
    }
    public static Matriz bidiagonalInferior(BigComplex[] diagonal, BigComplex[] inferior) {
        return bidiagonalSuperior(diagonal, inferior).traspuesta();
    }
    public static Matriz tridiagonal(BigComplex[] inferior, BigComplex[] diagonal, BigComplex[] superior) {
        if (diagonal == null || inferior == null || superior == null
                || inferior.length != diagonal.length - 1 || superior.length != diagonal.length - 1)
            throw new DimensionMatrizException("Bandas incompatibles");
        Matriz m = bidiagonalSuperior(diagonal, superior);
        for (int i = 1; i < diagonal.length; i++) m = m.conElemento(i, i - 1, inferior[i - 1]);
        return m;
    }
    public static Matriz banda(int filas, int columnas, int inferior, int superior, BigComplex valor) {
        if (inferior < 0 || superior < 0) throw new DimensionMatrizException("Ancho de banda negativo");
        BigComplex[][] d = datosDimensiones(filas, columnas, BigComplex.ZERO, false);
        for (int i = 0; i < filas; i++) for (int j = 0; j < columnas; j++)
            if (i - j <= inferior && j - i <= superior) d[i][j] = Objects.requireNonNull(valor);
        return new Matriz(d);
    }
    public static Matriz hessenbergSuperior(BigComplex[][] datos) {
        Matriz m = new Matriz(datos);
        if (!m.esHessenbergSuperior()) throw new DimensionMatrizException("No es Hessenberg superior");
        return m;
    }
    public static Matriz desdeMatrizAdyacencia(BigComplex[][] datos) {
        Matriz m = new Matriz(datos);
        if (!m.esCuadrada() || !m.esNoNegativa()) throw new DimensionMatrizException("Adyacencia inválida");
        return m;
    }
    public static Matriz matrizDeIncidencia(int vertices, int[][] aristas, boolean dirigida) {
        if (vertices <= 0 || aristas == null) throw new DimensionMatrizException("Grafo inválido");
        BigComplex[][] d = datosDimensiones(vertices, Math.max(1, aristas.length), BigComplex.ZERO, false);
        for (int e = 0; e < aristas.length; e++) {
            if (aristas[e] == null || aristas[e].length != 2 || aristas[e][0] < 0 || aristas[e][0] >= vertices
                    || aristas[e][1] < 0 || aristas[e][1] >= vertices) throw new DimensionMatrizException("Arista inválida");
            d[aristas[e][0]][e] = dirigida ? Complejos.MINUS_ONE : BigComplex.ONE;
            d[aristas[e][1]][e] = BigComplex.ONE;
        }
        return new Matriz(d);
    }
    public static Matriz matrizCompanera(BigComplex... coeficientesMonicos) {
        if (coeficientesMonicos == null || coeficientesMonicos.length < 2
                || !Complejos.aproximadamenteIguales(coeficientesMonicos[0], BigComplex.ONE,
                Complejos.TOLERANCIA_DEFECTO, Complejos.CONTEXTO_DEFECTO))
            throw new DimensionMatrizException("Se requieren coeficientes mónicos en orden descendente");
        int n = coeficientesMonicos.length - 1;
        BigComplex[][] d = datosDimensiones(n, n, BigComplex.ZERO, false);
        for (int i = 1; i < n; i++) d[i][i - 1] = BigComplex.ONE;
        for (int i = 0; i < n; i++) d[i][n - 1] = coeficientesMonicos[n - i].negate();
        return new Matriz(d);
    }
    public static Matriz bloqueJordan(BigComplex lambda, int n) {
        BigComplex[][] d = datosDimensiones(n, n, BigComplex.ZERO, false);
        for (int i = 0; i < n; i++) {
            d[i][i] = Objects.requireNonNull(lambda);
            if (i + 1 < n) d[i][i + 1] = BigComplex.ONE;
        }
        return new Matriz(d);
    }

    private static BigComplex[][] datosDimensiones(int f, int c, BigComplex valor, boolean diagonal) {
        if (f <= 0 || c <= 0) throw new DimensionMatrizException("Dimensiones positivas requeridas");
        Objects.requireNonNull(valor);
        BigComplex[][] d = new BigComplex[f][c];
        for (int i = 0; i < f; i++) for (int j = 0; j < c; j++)
            d[i][j] = diagonal ? (i == j ? BigComplex.ONE : BigComplex.ZERO) : valor;
        return d;
    }
    private static BigComplex[][] convertir(double[][] a) {
        Objects.requireNonNull(a); BigComplex[][] d = new BigComplex[a.length][];
        for (int i = 0; i < a.length; i++) {
            if (a[i] == null) throw new DimensionMatrizException("Fila null");
            d[i] = new BigComplex[a[i].length];
            for (int j = 0; j < a[i].length; j++) d[i][j] = Complejos.de(a[i][j]);
        }
        return d;
    }
    private static BigComplex[][] convertir(Double[][] a) {
        Objects.requireNonNull(a); BigComplex[][] d = new BigComplex[a.length][];
        for (int i = 0; i < a.length; i++) {
            if (a[i] == null) throw new DimensionMatrizException("Fila null");
            d[i] = new BigComplex[a[i].length];
            for (int j = 0; j < a[i].length; j++) d[i][j] = Complejos.de(Objects.requireNonNull(a[i][j]));
        }
        return d;
    }
    private static BigComplex[][] convertir(BigDecimal[][] a) {
        Objects.requireNonNull(a); BigComplex[][] d = new BigComplex[a.length][];
        for (int i = 0; i < a.length; i++) {
            if (a[i] == null) throw new DimensionMatrizException("Fila null");
            d[i] = new BigComplex[a[i].length];
            for (int j = 0; j < a[i].length; j++) d[i][j] = BigComplex.valueOf(Objects.requireNonNull(a[i][j]));
        }
        return d;
    }
    private static BigComplex[][] convertir(String[][] a) {
        Objects.requireNonNull(a); BigComplex[][] d = new BigComplex[a.length][];
        for (int i = 0; i < a.length; i++) {
            if (a[i] == null) throw new DimensionMatrizException("Fila null");
            d[i] = new BigComplex[a[i].length];
            for (int j = 0; j < a[i].length; j++) d[i][j] = Complejos.de(a[i][j]);
        }
        return d;
    }
    private static BigComplex[][] validar(BigComplex[][] a, MathContext mc) {
        Objects.requireNonNull(a);
        if (a.length == 0 || a[0] == null || a[0].length == 0) throw new DimensionMatrizException("Matriz vacía");
        int c = a[0].length; BigComplex[][] d = new BigComplex[a.length][c];
        for (int i = 0; i < a.length; i++) {
            if (a[i] == null || a[i].length != c) throw new DimensionMatrizException("Matriz irregular");
            for (int j = 0; j < c; j++) d[i][j] = Objects.requireNonNull(a[i][j]).round(mc);
        }
        return d;
    }
    private static BigComplex[][] copiar(MatrizCompleja m) {
        Objects.requireNonNull(m); BigComplex[][] d = new BigComplex[m.filas()][m.columnas()];
        for (int i = 0; i < m.filas(); i++) for (int j = 0; j < m.columnas(); j++) d[i][j] = m.get(i, j);
        return d;
    }
    private static BigComplex[][] convertirDouble(ArrayList<ArrayList<Double>> v) {
        Objects.requireNonNull(v);
        if (v.isEmpty() || v.getFirst() == null || v.getFirst().isEmpty()) throw new DimensionMatrizException("Lista vacía");
        BigComplex[][] d = new BigComplex[v.size()][v.getFirst().size()];
        for (int i = 0; i < v.size(); i++) {
            if (v.get(i) == null || v.get(i).size() != d[0].length) throw new DimensionMatrizException("Lista irregular");
            for (int j = 0; j < d[i].length; j++) d[i][j] = Complejos.de(Objects.requireNonNull(v.get(i).get(j)));
        }
        return d;
    }
    private static BigComplex[][] convertirNumber(List<? extends List<? extends Number>> v) {
        Objects.requireNonNull(v);
        if (v.isEmpty() || v.getFirst() == null || v.getFirst().isEmpty()) throw new DimensionMatrizException("Lista vacía");
        BigComplex[][] d = new BigComplex[v.size()][v.getFirst().size()];
        for (int i = 0; i < v.size(); i++) {
            if (v.get(i) == null || v.get(i).size() != d[0].length) throw new DimensionMatrizException("Lista irregular");
            for (int j = 0; j < d[i].length; j++) d[i][j] = Complejos.de(Objects.requireNonNull(v.get(i).get(j)));
        }
        return d;
    }
}
