FROM maven:3.9.11-eclipse-temurin-21 AS build
# CVE-2026-114938: temporarily disable the vulnerable Base64 intrinsic.
ENV JAVA_TOOL_OPTIONS="-XX:-UseBASE64Intrinsics"
WORKDIR /app
COPY pom.xml .
RUN mvn -q dependency:go-offline
COPY src src
RUN mvn -q clean test package dependency:copy-dependencies -DoutputDirectory=target/lib

FROM eclipse-temurin:21-jre
ENV JAVA_TOOL_OPTIONS="-XX:-UseBASE64Intrinsics"
WORKDIR /app
COPY --from=build /app/target/biblioteca-matrices-1.0.0.jar app.jar
COPY --from=build /app/target/lib lib
ENTRYPOINT ["java", "-cp", "app.jar:lib/*", "BibliotecaMatricesApp"]
