# Visualizing Model Performance

## Models are easier to communicate with interesting and informative graphics. In this repository, I share a few data visualizations that have helped me get my message across.

### Background
Imagine you are an ecologist studying the Sulu woodpecker (*Yungipicus ramsayi*). Populations are decreasing according to the [IUCN](https://www.iucnredlist.org/fr/search/list?taxonomies=22673150&searchType=species). It is up to you to learn as much as you can about this species to determine how you can help restore native populations. 

You decide that predictive modeling could help you identify the most important habitat characteristics for this species and identify geographical areas where it could thrive. First, you collect data in three types of ecosystems where you know Sulu woodpecker occur: a lowland forest, a mangrove forest, and cultivated cropland. Next, you gather predictor variables that best represent this species' ecological needs. These include climate variables (e.g., solar radiation, temperature, precipitation), biotic variables (e.g., the presence of specific trees, predator abundance, food source abundance), and human variables (e.g., nearby road density, level of human disturbance). Finally, you are ready to model your data across these ecosystems to learn about Sulu woodpecker needs and predict what areas could contain the best habitat. You don't know which variables are the most important biologically, so you decide to try different combinations of climate, biotic, and human predictors.

Now that you have your models, you're not sure which ones to use! Luckily, you collected plenty of data and you're able to build multiple models for each ecosystem and combination of predictors. That way, you can test the average predictive performance of each combination of variables in each ecosystem. You select three metrics to assess predictive performance: 1) Precision, which is the proporition of positive predictions that are actually positives; 2) Recall, which is the proportion of actually positives that are correctly predicted as positives; and 3) F1 Score, which balance Precision and Recall into one metric.

Looking at tables of Precision, Recall, and F1 Score is confusing. It's time to visualize your performance metrics to learn more about your models' performance!

### Explanation
Each subfolder demonstrates a different method for visualizing model performance. The folders contain R code for making the data visualizations, associated materials like custom icons, and example plots to help you check that your code is doing what you expect it to do. 

Best of luck to you and to the Sulu woodpecker!

## Note: This is fake data I created for the purpose of this project. I am not a Sulu woodpecker expert. The data are fabricated and should not be used as evidence of any particular belief about Sulu woodpeckers. If you would like examples of data visualization in an actual biodiversity project, check out my [OSF repository](https://osf.io/y6meq/) on data-centric species distribution modeling of invasive European frog-bit in North America.



