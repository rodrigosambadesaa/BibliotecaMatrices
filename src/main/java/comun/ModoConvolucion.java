package comun;

/** Tamaño de salida para convoluciones y correlaciones tridimensionales. */
public enum ModoConvolucion {
    /** Todas las superposiciones, incluyendo las parciales con relleno cero. */
    COMPLETA,
    /** Salida con las mismas dimensiones que el primer operando. */
    MISMA,
    /** Solo posiciones donde el núcleo queda completamente dentro del tensor. */
    VALIDA
}
