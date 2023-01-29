-- MySQL dump 10.13  Distrib 8.0.28, for Win64 (x86_64)
--
-- Host: localhost    Database: frituurdb
-- ------------------------------------------------------
-- Server version	8.0.28

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!50503 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `reviews`
--

DROP TABLE IF EXISTS `reviews`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `reviews` (
  `id_review` int NOT NULL AUTO_INCREMENT,
  `id_author` int NOT NULL,
  `score` double(10,3) DEFAULT NULL,
  `id_frituur` int DEFAULT NULL,
  PRIMARY KEY (`id_review`),
  KEY `id_author` (`id_author`),
  KEY `id_frituur` (`id_frituur`),
  CONSTRAINT `reviews_ibfk_1` FOREIGN KEY (`id_author`) REFERENCES `review_authors` (`id_author`),
  CONSTRAINT `reviews_ibfk_2` FOREIGN KEY (`id_frituur`) REFERENCES `frituur` (`id_frituur`)
) ENGINE=InnoDB AUTO_INCREMENT=30 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_0900_ai_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `reviews`
--

LOCK TABLES `reviews` WRITE;
/*!40000 ALTER TABLE `reviews` DISABLE KEYS */;
INSERT INTO `reviews` VALUES (1,1,6.000,76),(2,2,5.000,76),(3,4,6.500,76),(4,5,7.123,76),(5,1,7.000,20),(6,2,7.500,20),(7,1,9.500,82),(8,2,9.500,82),(9,1,6.000,55),(10,2,6.000,55),(11,1,6.500,38),(12,2,6.500,38),(13,1,7.500,83),(14,2,7.500,83),(15,1,7.700,78),(16,2,7.500,78),(17,3,8.200,78),(18,4,8.000,78),(19,5,7.500,78),(20,1,8.050,18),(21,2,8.000,18),(22,3,9.000,18),(23,4,7.500,18),(24,5,8.250,18),(25,1,8.200,81),(26,2,8.500,81),(27,3,8.000,81),(28,4,7.750,81),(29,5,7.950,81);
/*!40000 ALTER TABLE `reviews` ENABLE KEYS */;
UNLOCK TABLES;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2022-05-01 10:45:58
