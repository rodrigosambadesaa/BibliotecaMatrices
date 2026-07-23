import ch.obermuhlner.math.big.BigComplex;
import comun.AbstractMatriz;
import comun.Complejos;
import comun.MatrizCompleja;
import comun.resultados.ResultadoAutovalor;
import comun.resultados.ResultadoDiagonalizacion;
import comun.resultados.ResultadoIterativo;
import comun.resultados.ResultadoLU;
import comun.resultados.ResultadoQR;
import comun.resultados.ResultadoSVD;
import comun.resultados.ResultadoSistema;
import otros.MatrizBooleana;
import otros.MatrizModular;

import java.lang.reflect.Array;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.function.BinaryOperator;
import java.util.function.Function;
import java.util.function.UnaryOperator;

/**
 * Demostración autocontenida. Además de las comprobaciones algebraicas
 * explícitas, descubre y ejecuta individualmente cada método público de las dos
 * matrices; muestra nombre, entradas, salida y validación, y falla si alguno no
 * pudo invocarse.
 */
public final class BibliotecaMatricesApp {
    private static final BigDecimal TOL = new BigDecimal("1E-18");
    private static int descubiertos;
    private static int invocados;
    private static int superados;
    private static final Map<String, String> fallidos = new LinkedHashMap<>();

    private BibliotecaMatricesApp() {
    }

    public static void main(String[] args) {
        System.out.println("BIBLIOTECA MATRICIAL BIGCOMPLEX");
        probarConstructoresArray();
        probarConstructoresLista();
        probarPropiedadesAlgebraicas();
        probarFactorizaciones();
        probarMatricesBooleanasYModulares();
        auditarApi(array.Matriz.class);
        auditarApi(lista.Matriz.class);
        auditarResultados();
        int noProbados = descubiertos - invocados;
        System.out.printf("%nRESUMEN: descubiertos=%d, invocados=%d, superados=%d, fallidos=%d, no probados=%d%n",
                descubiertos, invocados, superados, fallidos.size(), noProbados);
        fallidos.forEach((metodo, causa) -> System.out.println("  FALLO " + metodo + ": " + causa));
        if (!fallidos.isEmpty() || noProbados != 0) {
            throw new IllegalStateException("La demostración exhaustiva no fue superada");
        }
    }

    private static void probarConstructoresArray() {
        seccion("Constructores array");
        probarConstructores(array.Matriz.class);
    }

    private static void probarConstructoresLista() {
        seccion("Constructores lista");
        probarConstructores(lista.Matriz.class);
    }

    private static void probarConstructores(Class<?> tipo) {
        for (Constructor<?> c : tipo.getConstructors()) {
            descubiertos++;
            String id = tipo.getName() + constructorFirma(c);
            try {
                Object[] entrada = argumentos(c.getParameterTypes(), c.getParameters(), c.getName(), null);
                System.out.println(id + " entrada=" + Arrays.deepToString(entrada));
                Object salida = c.newInstance(entrada);
                comprobar(salida instanceof MatrizCompleja m && m.filas() > 0 && m.columnas() > 0,
                        "El constructor no produjo una matriz válida");
                invocados++; superados++;
                System.out.println("  resultado=" + resumen(salida) + " CORRECTO");
            } catch (ReflectiveOperationException | RuntimeException e) {
                invocados++; fallidos.put(id, causa(e));
            }
        }
    }

    private static void probarPropiedadesAlgebraicas() {
        seccion("Propiedades algebraicas cruzadas");
        array.Matriz a = new array.Matriz(new double[][]{{4, 1}, {1, 3}});
        array.Matriz b = new array.Matriz(new double[][]{{1, 2}, {0, 1}});
        comprobar(a.suma(array.Matriz.cero(2, 2)).equals(a), "A+0=A");
        comprobar(a.resta(a).esNula(), "A-A=0");
        comprobar(a.producto(b.suma(a)).aproximadamenteIgual(a.producto(b).suma(a.producto(a)), TOL),
                "Distributividad");
        comprobar(a.producto(b).traspuesta().aproximadamenteIgual(b.traspuesta().producto(a.traspuesta()), TOL),
                "Traspuesta del producto");
        comprobar(a.verificarInversa(a.inversa(), TOL), "Inversa");
        lista.Matriz l = new lista.Matriz(a);
        comprobar(a.aproximadamenteIgual(l, TOL), "Equivalencia array/lista");
        comprobar(l.suma(l) instanceof lista.Matriz, "La implementación lista conserva su tipo");
        System.out.println("Identidades básicas y equivalencia de almacenamientos: CORRECTO");
    }

    private static void probarFactorizaciones() {
        seccion("Factorizaciones y residuos");
        array.Matriz a = new array.Matriz(new double[][]{{4, 1}, {1, 3}});
        comprobar(a.verificarLU(a.descomposicionLUP(), TOL), "PA=LU");
        comprobar(a.verificarQR(a.descomposicionQRGramSchmidtModificado(), TOL), "A=QR");
        comprobar(a.verificarCholesky(a.descomposicionCholesky(), TOL), "A=LL*");
        comprobar(a.resolverJacobi(array.Matriz.columna(BigComplex.ONE, BigComplex.ONE), TOL, 200).convergio(),
                "Jacobi");
        System.out.println("LU, QR, Cholesky e iteración: CORRECTO");
    }

    private static void probarMatricesBooleanasYModulares() {
        seccion("Otros cuerpos");
        MatrizBooleana b = new MatrizBooleana(new boolean[][]{{false, true}, {false, false}});
        comprobar(b.clausuraTransitiva().get(0, 1), "Clausura booleana");
        MatrizModular m = new MatrizModular(new int[][]{{1, 2}, {3, 4}}, 5);
        comprobar(m.producto(m.inversa()).equals(new MatrizModular(new int[][]{{1, 0}, {0, 1}}, 5)),
                "Inversa modular");
        System.out.println("Semianillo booleano y cuerpo primo: CORRECTO");
    }

    private static void auditarApi(Class<?> tipo) {
        seccion("Auditoría reflectiva " + tipo.getName());
        Method[] metodos = Arrays.stream(tipo.getMethods())
                .filter(BibliotecaMatricesApp::metodoDeBiblioteca)
                .sorted(Comparator.comparing(Method::getName).thenComparing(Method::toGenericString))
                .toArray(Method[]::new);
        for (Method metodo : metodos) {
            descubiertos++;
            String id = tipo.getName() + "#" + metodo.getName() + firma(metodo);
            Object receptor = Modifier.isStatic(metodo.getModifiers()) ? null : receptorPara(tipo, metodo.getName());
            try {
                Object[] entrada = argumentos(metodo.getParameterTypes(), metodo.getParameters(), metodo.getName(), receptor);
                System.out.println(id + " entrada=" + Arrays.deepToString(entrada));
                Object salida = metodo.invoke(receptor, entrada);
                comprobar(metodo.getReturnType() == void.class || salida != null,
                        "Un método público devolvió null");
                validarResultado(metodo, receptor, salida);
                invocados++; superados++;
                System.out.println("  resultado=" + resumen(salida) + " CORRECTO");
            } catch (ReflectiveOperationException | RuntimeException e) {
                invocados++; fallidos.put(id, causa(e));
            }
        }
    }

    private static boolean metodoDeBiblioteca(Method m) {
        Package p = m.getDeclaringClass().getPackage();
        return Modifier.isPublic(m.getModifiers()) && p != null
                && (p.getName().equals("array") || p.getName().equals("lista") || p.getName().equals("comun"));
    }

    private static Object receptorPara(Class<?> tipo, String nombre) {
        if (nombre.contains("Transicion") || nombre.equals("distribuirTrasPasos")
                || nombre.equals("distribucionEstacionaria") || nombre.equals("esCadenaRegular")) {
            return tipo == array.Matriz.class
                    ? new array.Matriz(new double[][]{{0.5, 0.5}, {0.25, 0.75}})
                    : new lista.Matriz(new double[][]{{0.5, 0.5}, {0.25, 0.75}});
        }
        if (nombre.contains("Diagonal") && nombre.startsWith("inversa")) {
            return tipo == array.Matriz.class ? array.Matriz.diagonal(Complejos.de(2), Complejos.de(3))
                    : lista.Matriz.diagonal(Complejos.de(2), Complejos.de(3));
        }
        if (nombre.equals("inversaTriangular")) {
            return tipo == array.Matriz.class ? array.Matriz.diagonal(Complejos.de(2), Complejos.de(3))
                    : lista.Matriz.diagonal(Complejos.de(2), Complejos.de(3));
        }
        if (nombre.equals("determinante3x3Sarrus")) {
            return tipo == array.Matriz.class ? array.Matriz.identidad(3) : lista.Matriz.identidad(3);
        }
        if (nombre.equals("indiceDeNilpotencia")) {
            return tipo == array.Matriz.class ? new array.Matriz(new double[][]{{0, 1}, {0, 0}})
                    : new lista.Matriz(new double[][]{{0, 1}, {0, 0}});
        }
        if (nombre.equals("periodo")) {
            return tipo == array.Matriz.class ? array.Matriz.identidad(2) : lista.Matriz.identidad(2);
        }
        if (nombre.equals("desvectorizar")) {
            MatrizCompleja a = matrizSpd(tipo);
            return tipo == array.Matriz.class ? a.comoMatrizArray().vectorizarPorFilas()
                    : a.comoMatrizLista().vectorizarPorFilas();
        }
        return matrizSpd(tipo);
    }

    private static MatrizCompleja matrizSpd(Class<?> tipo) {
        return tipo == array.Matriz.class
                ? new array.Matriz(new double[][]{{4, 1}, {1, 3}})
                : new lista.Matriz(new double[][]{{4, 1}, {1, 3}});
    }

    private static Object[] argumentos(Class<?>[] tipos, java.lang.reflect.Parameter[] parametros,
                                       String nombre, Object receptor) {
        Object[] args = new Object[tipos.length];
        for (int i = 0; i < tipos.length; i++) args[i] = argumento(tipos[i], parametros[i].getName(), nombre, i, receptor);
        return args;
    }

    private static Object argumento(Class<?> t, String parametro, String metodo, int posicion, Object receptor) {
        Class<?> concreto = receptor instanceof lista.Matriz || metodo.contains("lista.Matriz")
                ? lista.Matriz.class : array.Matriz.class;
        if (t == int.class || t == Integer.class) return entero(parametro, metodo, posicion);
        if (t == long.class || t == Long.class) return 2L;
        if (t == double.class || t == Double.class) return 2.0;
        if (t == boolean.class || t == Boolean.class) return false;
        if (t == char.class || t == Character.class) return 'z';
        if (t == BigDecimal.class) {
            if (parametro.toLowerCase().contains("tol") || parametro.equals("paso")) return TOL;
            if (parametro.toLowerCase().contains("omega")) return BigDecimal.ONE;
            if (parametro.toLowerCase().contains("densidad")) return new BigDecimal("0.2");
            return BigDecimal.ONE;
        }
        if (t == MathContext.class) return MathContext.DECIMAL128;
        if (t == BigComplex.class) return metodo.contains("lambda") ? Complejos.de(3) : BigComplex.ONE;
        if (t == BigComplex[].class) {
            if (metodo.startsWith("bidiagonal")) return posicion == 0
                    ? new BigComplex[]{Complejos.de(2), Complejos.de(3)}
                    : new BigComplex[]{BigComplex.ONE};
            if (metodo.equals("tridiagonal")) return posicion == 1
                    ? new BigComplex[]{Complejos.de(2), Complejos.de(3)}
                    : new BigComplex[]{BigComplex.ONE};
            if (metodo.equals("hankel")) return posicion == 0
                    ? new BigComplex[]{BigComplex.ONE, Complejos.de(2)}
                    : new BigComplex[]{Complejos.de(2), Complejos.de(3)};
            return new BigComplex[]{BigComplex.ONE, BigComplex.ONE};
        }
        if (t == BigComplex[][].class) {
            if (metodo.equals("triangularSuperior")) return new BigComplex[][]{{Complejos.de(4), BigComplex.ONE},
                    {BigComplex.ZERO, Complejos.de(3)}};
            if (metodo.equals("triangularInferior")) return new BigComplex[][]{{Complejos.de(4), BigComplex.ZERO},
                    {BigComplex.ONE, Complejos.de(3)}};
            return new BigComplex[][]{{Complejos.de(4), BigComplex.ONE}, {BigComplex.ONE, Complejos.de(3)}};
        }
        if (t == BigDecimal[][].class) return new BigDecimal[][]{{BigDecimal.valueOf(4), BigDecimal.ONE},
                {BigDecimal.ONE, BigDecimal.valueOf(3)}};
        if (t == double[][].class) return new double[][]{{4, 1}, {1, 3}};
        if (t == Double[][].class) return new Double[][]{{4.0, 1.0}, {1.0, 3.0}};
        if (t == String[][].class) return new String[][]{{"4", "1"}, {"1", "3"}};
        if (t == int[].class) {
            if (metodo.equals("matrizConfusion")) return new int[]{0, 1};
            if (parametro.toLowerCase().contains("taman")) return new int[]{1, 1};
            if (metodo.contains("extraer") || metodo.contains("permutar")) return new int[]{0, 1};
            if (metodo.contains("Permutacion")) return new int[]{1, 0};
            return new int[]{1};
        }
        if (MatrizCompleja.class.isAssignableFrom(t) || AbstractMatriz.class.isAssignableFrom(t)) {
            boolean iterativo = metodo.equals("GMRES") || metodo.equals("BiCG") || metodo.equals("BiCGSTAB")
                    || metodo.equals("gradienteConjugado") || metodo.equals("resolverJacobi")
                    || metodo.equals("resolverGaussSeidel") || metodo.equals("resolverSOR");
            boolean sistemaVector = metodo.equals("resolver") || metodo.equals("resolverPorGauss")
                    || metodo.equals("resolverPorGaussJordan") || metodo.equals("resolverPorLU")
                    || metodo.equals("resolverPorQR") || metodo.equals("resolverPorCholesky")
                    || metodo.equals("resolverPorCramer") || metodo.equals("minimosCuadrados")
                    || metodo.equals("minimosCuadradosQR") || metodo.equals("minimosCuadradosSVD")
                    || metodo.equals("ecuacionesNormales") || metodo.equals("resolverConPseudoinversa")
                    || metodo.equals("regresionMinimosCuadrados") || metodo.equals("refinamientoIterativo");
            MatrizCompleja m = posicion == 0 && (iterativo || sistemaVector)
                    ? vector(concreto) : matrizSpd(concreto);
            return t == lista.Matriz.class ? m.comoMatrizLista()
                    : t == array.Matriz.class ? m.comoMatrizArray() : m;
        }
        if (t == ArrayList.class) {
            if (metodo.equals("desdeBigComplex")) {
                ArrayList<ArrayList<BigComplex>> v = new ArrayList<>();
                v.add(new ArrayList<>(List.of(Complejos.de(4), BigComplex.ONE)));
                v.add(new ArrayList<>(List.of(BigComplex.ONE, Complejos.de(3))));
                return v;
            }
            ArrayList<ArrayList<Double>> v = new ArrayList<>();
            v.add(new ArrayList<>(List.of(4.0, 1.0))); v.add(new ArrayList<>(List.of(1.0, 3.0))); return v;
        }
        if (t == List.class) {
            if (metodo.equals("lista.Matriz")) return List.of(List.of(4, 1), List.of(1, 3));
            if (metodo.equals("aplicarOperacionesElementales")) return List.of();
            return List.of(matrizSpd(concreto));
        }
        if (t == UnaryOperator.class) return (UnaryOperator<BigComplex>) z -> z;
        if (t == BinaryOperator.class) return (BinaryOperator<BigComplex>) BigComplex::add;
        if (t == Function.class) return (Function<BigDecimal, MatrizCompleja>) x ->
                concreto == lista.Matriz.class
                        ? matrizSpd(concreto).comoMatrizLista().productoEscalar(x)
                        : matrizSpd(concreto).comoMatrizArray().productoEscalar(x);
        if (t == MatrizCompleja[][].class) {
            MatrizCompleja uno = concreto == lista.Matriz.class ? lista.Matriz.identidad(1) : array.Matriz.identidad(1);
            return new MatrizCompleja[][]{{uno, uno}, {uno, uno}};
        }
        if (t == ResultadoLU.class) return matrizSpd(concreto).comoMatrizArray().descomposicionLUP();
        if (t == ResultadoQR.class) return matrizSpd(concreto).comoMatrizArray().descomposicionQRGramSchmidtModificado();
        if (t == ResultadoSVD.class) return array.Matriz.diagonal(Complejos.de(2), Complejos.de(3)).descomposicionSVD();
        if (t == ResultadoDiagonalizacion.class) return array.Matriz.diagonal(Complejos.de(2), Complejos.de(3)).diagonalizacion();
        if (t.isArray()) return Array.newInstance(t.componentType(), 0);
        if (t.isAssignableFrom(concreto)) return matrizSpd(concreto);
        throw new IllegalArgumentException("No hay dato de demostración para " + t.getTypeName());
    }

    private static int entero(String p, String metodo, int posicion) {
        String n = p.toLowerCase();
        if (n.contains("max") || n.contains("iteracion")) return 40;
        if (n.equals("filas") || n.equals("columnas") || n.equals("f") || n.equals("c")) return 2;
        if (n.contains("vertices")) return 2;
        if (n.contains("intervalos")) return 20;
        if (n.contains("rango")) return 1;
        if (n.contains("clases")) return 2;
        if (n.contains("pasos")) return 2;
        if (metodo.equals("potenciaNegativa")) return -2;
        if (n.contains("exponente") || n.contains("longitud") || n.equals("k")) return 2;
        if (n.contains("nuevasfilas") || n.contains("nuevascolumnas")) return 2;
        if (n.contains("filahasta") || n.contains("columnahasta")) return 1;
        if (n.contains("filadesde") || n.contains("columnadesde")) return 0;
        if (metodo.contains("intercambiar")) return posicion;
        if (n.equals("j") || n.contains("destino") || n.contains("columna")) return 1;
        if (n.equals("n")) return 2;
        if (n.contains("corte")) return 1;
        if (n.contains("escala")) return 3;
        return 0;
    }

    private static MatrizCompleja vector(Class<?> concreto) {
        return concreto == lista.Matriz.class
                ? lista.Matriz.columna(BigComplex.ONE, BigComplex.ONE)
                : array.Matriz.columna(BigComplex.ONE, BigComplex.ONE);
    }

    private static void validarResultado(Method metodo, Object receptor, Object salida) {
        if (salida instanceof MatrizCompleja m) comprobar(m.filas() > 0 && m.columnas() > 0, "Dimensión de salida inválida");
        if (salida instanceof BigDecimal d) comprobar(d.signum() >= 0 || metodo.getName().contains("angulo"),
                "Magnitud inesperadamente negativa");
        if (salida instanceof ResultadoIterativo<?> r) comprobar(r.iteraciones() >= 0, "Iteraciones negativas");
        if (salida instanceof ResultadoAutovalor<?> r) comprobar(r.residuo().signum() >= 0, "Residuo negativo");
        if (receptor instanceof MatrizCompleja original && salida instanceof MatrizCompleja resultado
                && metodo.getName().equals("traspuesta")) {
            comprobar(resultado.filas() == original.columnas(), "Dimensiones de transpuesta incorrectas");
        }
    }

    private static void auditarResultados() {
        seccion("Clases públicas de resultado");
        array.Matriz a = new array.Matriz(new double[][]{{4, 1}, {1, 3}});
        List<Object> resultados = List.of(
                a.formaEscalonadaConDetalles(), a.descomposicionLUP(),
                a.descomposicionQRGramSchmidtModificado(),
                a.metodoDeLaPotencia(TOL, 20), a.discutirSistema(array.Matriz.columna(BigComplex.ONE, BigComplex.ONE)));
        for (Object resultado : resultados) {
            for (Method m : resultado.getClass().getMethods()) {
                if (m.getParameterCount() == 0 && m.getDeclaringClass().getPackageName().equals("comun.resultados")) {
                    descubiertos++;
                    try {
                        Object valor = m.invoke(resultado);
                        comprobar(valor != null, "Accessor de resultado devolvió null");
                        invocados++; superados++;
                    } catch (ReflectiveOperationException | RuntimeException e) {
                        invocados++; fallidos.put(resultado.getClass().getName() + "#" + m.getName(), causa(e));
                    }
                }
            }
        }
    }

    private static void comprobar(boolean condicion, String mensaje) {
        if (!condicion) throw new IllegalStateException(mensaje);
    }
    private static String constructorFirma(Constructor<?> c) {
        return Arrays.toString(c.getParameterTypes());
    }
    private static String firma(Method m) { return Arrays.toString(m.getParameterTypes()); }
    private static String resumen(Object o) {
        if (o == null) return "void";
        String s = Objects.toString(o).replace(System.lineSeparator(), " / ");
        return s.length() > 160 ? s.substring(0, 157) + "..." : s;
    }
    private static String causa(Throwable e) {
        Throwable t = e instanceof InvocationTargetException i && i.getCause() != null ? i.getCause() : e;
        return t.getClass().getSimpleName() + ": " + t.getMessage();
    }
    private static void seccion(String nombre) { System.out.println("\n== " + nombre + " =="); }
}
