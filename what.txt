# **Advanced Minecraft Seed Finder**  

## **Website:** [mc-seed-finder.replit.app](https://mc-seed-finder.replit.app/)  

### *(If website is down, you can use the seed finder by following the Installation)*

## **Overview**  
This tool automates finding Minecraft seeds based on specific biomes, structures, and terrain features. Supports Java and Bedrock editions (1.18 to 1.21.6)

## **Features**  
- Find any biomes including custom biomes (like islands, valleys, and elevated encircling terrain)
- Find structures with biome constraints, height ranges, and minimum counts  
- Locate clustered structures (structures next to each other) and biome combinations  
- Set biome size requirements for more precise world generation  

## **Usage**  
1. Define search parameters  
2. Select Minecraft version and edition  
3. Set search center using coordinates or spawn point
4. Set the seed requirments (such as at least 1 village and Cherry Grove within 500 block search radius from spawn for example)
6. Run the scan, might take a few minutes depending on the rarity of the seed requirements

## **Installation**  

### **Clone the Repository**  
```bash
git clone https://github.com/Dalle2fan/cubiomes.git
cd cubiomes
```

```bash
brew install cmake
```

### **Run Application**  
```bash
./build.sh
``` 
