## Libraries -----------------------------------


.libPaths(new = "C:/Users/Stefon Flego/Programs/R/R-4.4.3/library")

 
# install.packages("colorspace", dependencies = TRUE)


library(tidyverse)
library(foreach)
library(mgcv)
library(gganimate)
library(av)



## I Indiana English b_d --------------------------------------------------------


dir1 <- "Indiana English vowels read in b_d context/"


# Read in data
data <- paste0(dir1, "Formants for American English male speaker producing vowels in b_d context.csv") %>% 
  read.csv(header = TRUE)


# Take a peak at the data set
View(data)


# specify factor level order for 'Vowel' variable (lexical set labels)
Vowel_levels <- c("CHOICE", "DRESS", "FACE", "FLEECE", "FOOT", "GOAT", "GOOSE", "KIT",
  "LOT", "MOUTH", "MUTE", "PRICE", "STRUT", "THOUGHT", "TRAP")


# Pre-processing! The visualization technique works best when formant
# frequencies have been converted to your preferred (non-linear) transformation
# of Hz (e.g. log(Hz), Lobanov, bark, etc.) For simplicity here, I'm going to
# remove the first and last timepoint measurements taken for each vowel, and
# I'll transform Hz logarithmically.

data <- data %>% 
  group_by(Token_ID) %>% 
  filter(NormTime != min(NormTime) & NormTime != max(NormTime)) %>% 
  ungroup() %>% 
  mutate(F1_Hz = F1, F2_Hz = F2, F3_Hz = F3,
         F1 = log(F1_Hz), F2 = log(F2_Hz), F3 = log(F3_Hz)) %>% 
  mutate(Vowel = factor(Vowel, levels = Vowel_levels))


# look at raw formant measurements output from Praat
data %>% 
  ggplot(aes(x = IntervalTime*1000, group = Token_ID)) +
  geom_line(aes(y = F1)) +
  geom_line(aes(y = F2)) +
  geom_line(aes(y = F3)) +
  scale_y_continuous("Formant frequency (logHz)") +
  theme_bw(base_size = 12) +
  scale_x_continuous("Time (ms)") +
  facet_wrap(vars(Vowel, Word))



# Create smooth average trajectory for each vowel class using GAMs
get_smooth <- function(df, k) { # df is a group-split data frame, k is smoothing parameter
  
  foreach(Vwls = df, .combine = rbind) %do% {
              
              
    
              pred_F1 <-  gam(F1 ~ s(NormTime, k = k), data = Vwls)
              pred_F2 <-  gam(F2 ~ s(NormTime, k = k), data = Vwls)
              pred_F3 <-  gam(F3 ~ s(NormTime, k = k), data = Vwls)
              
              
              Vwls$Duration <- Vwls %>% 
                group_by(Vowel, Token_ID, Duration) %>% 
                summarize( .groups = "drop") %>% 
                group_by(Vowel) %>%  
                summarize(Duration = mean(Duration), .groups = "drop") %>% 
                pull(Duration)
              
              
              
              smooths <- Vwls %>% 
                group_by(Vowel, Word) %>% 
                reframe(
                  Duration = mean(unique(Duration)),
                  RealTime = seq( Duration*0, Duration*1, 0.001)
                ) %>% 
                mutate(NormTime = RealTime/Duration, .after = RealTime)
              
              
              smooths$F1 <- predict(pred_F1, newdata = smooths, type = "response")
              smooths$F2 <- predict(pred_F2, newdata = smooths, type = "response")
              smooths$F3 <- predict(pred_F3, newdata = smooths, type = "response")
              
              smooths
              
            }
  
}


smooth_data <- get_smooth(df = data %>% group_split(Vowel, Word), k = 5)



# Get velocity through F1 x F2 space at each time point for each vowel 

smooth_data <- smooth_data %>%   
  group_by(Vowel) %>% 
  mutate(F1_veloc = predict(smooth.spline(x = RealTime, y = F1), deriv = 1)$y,
         F2_veloc = predict(smooth.spline(x = RealTime, y = F2), deriv = 1)$y, 
         F3_veloc = predict(smooth.spline(x = RealTime, y = F3), deriv = 1)$y 
  ) %>% 
  ungroup() %>% 
  mutate(veloc = sqrt(F1_veloc^2 + F2_veloc^2)) # could add F3 to this if desired




# preliminaries for vowel Space Visualizations 


arrow_ends <-  smooth_data %>%
  arrange(desc(NormTime)) %>%
  group_by(Vowel) %>% slice(1) %>% ungroup()


# color palette used for vowel classes
{
  colors <- c("#D39200", # CHOICE
              "#00BA38", # DRESS
              "#93AA00", # FACE
              "#619CFF", # FLEECE 
              "#00C19F", # FOOT
              "#00ADFA", # GOAT
              "#B79F00", # GOOSE
              "#5EB300", # KIT
              "#FF699C", # LOT
              "#00B9E3", # MOUTH
              "#F8766D", # MUTE
              "#00BF74", # PRICE
              "#E88526", # STRUT
              "#DB72FB", # THOUGHT
              "#FF61C3"  # TRAP
  )
  
  }






## I.A plotting parameters that work well for a 15" vowel space image ------------------
## (like a conference poster) 


point_alpha <- 0.2
bgd_line_alpha <- point_alpha
arrow_b_size <- 10
max_point_size <- 20
min_point_size <- 2
bgd_line_width <- min_point_size + 1
point_stroke <- 0.35
breaks <- seq(0, 20, 0.25)

b_size = 50
theme_set(theme_bw(base_size = b_size))


# plot
p <- smooth_data %>%   
    arrange(NormTime, Vowel) %>% 
    ggplot(aes(x = F2, y = F1, group = Vowel, 
               color = Vowel, fill = Vowel)) +
    geom_path(linewidth = bgd_line_width, alpha = bgd_line_alpha) +
    geom_point(aes(size = veloc), color = "black", alpha = point_alpha, 
               shape = 21, stroke = point_stroke) +
    geom_point(data = arrow_ends, aes(size = veloc - arrow_b_size), shape = 24,
               color = "black", alpha = 0.75, stroke = point_stroke) +
    scale_x_reverse("F2 (logHz)", breaks = breaks) +
    scale_y_reverse("F1 (logHz)", breaks = breaks) +
    scale_fill_manual(values = colors 
                      # guide = "none"
                      ) +
    guides(fill = guide_legend(override.aes = list(size = 5))) +
    scale_color_manual(values = colors 
                       # guide = "none"
                       ) +
    scale_size_continuous(range = c(max_point_size, min_point_size), limits = c(-arrow_b_size, NA), guide = "none")
  
# save
  ggsave(p,
    filename = "American_English_vowels_b_d- large.png",
    path = dir1,
    dpi = "retina", width = unit(20, "inches"), height = unit(15, "inches"))
  
  

## I.B plotting parameters that work well for a large faceted image ------------------
## (each panel about 6")
  
  point_alpha <- 0.2
  bgd_line_alpha <- point_alpha
  arrow_b_size <- 10
  max_point_size <- 7
  min_point_size <- 1.5
  bgd_line_width <- min_point_size + 1
  point_stroke <- 0.25
  breaks <- seq(0, 20, 0.5)
  
  b_size = 50
  theme_set(theme_bw(base_size = b_size))
  
  
  # plot
  p <- smooth_data %>%   
    arrange(NormTime, Vowel) %>% 
    ggplot(aes(x = F2, y = F1, 
               color = Vowel, fill = Vowel)) +
    
    geom_point(data = smooth_data %>% select(-Vowel), aes(size = veloc), 
               shape = 21, fill = "white", color = "black", alpha = 0.05, stroke = point_stroke) +
    geom_point(data = arrow_ends %>% select(-Vowel), 
               aes(size = veloc - arrow_b_size),
               shape = 24, fill = "white", color = "black", alpha = 0.1, stroke = point_stroke) +
    
    geom_path(linewidth = bgd_line_width, alpha = bgd_line_alpha) +
    geom_point(aes(size = veloc), color = "black", 
               alpha = point_alpha, shape = 21, stroke = point_stroke) +
    geom_point(data = arrow_ends, aes(size = veloc - arrow_b_size), 
               shape = 24, color = "black", alpha = 0.75, stroke = point_stroke) +
    scale_x_reverse("F2 (logHz)", breaks = breaks) +
    scale_y_reverse("F1 (logHz)", breaks = breaks) +
    scale_fill_manual(values = colors, 
                      guide = "none"
    ) +
    # guides(fill = guide_legend(override.aes = list(size = 5))) +
    scale_color_manual(values = colors, 
                       guide = "none"
    ) +
    scale_size_continuous(range = c(max_point_size, min_point_size), limits = c(-arrow_b_size, NA), guide = "none") +
    facet_wrap(vars(Vowel), nrow = 3)
  
  # save
  ggsave(p,
         filename = "American_English_vowels_b_d- large faceted.png",
         path = dir1,
         dpi = "retina", width = unit(30, "inches"), height = unit(18, "inches"))
  
  
    
## I.C plotting parameters that work well for a 6" vowel space image ------------------
## (like in a paper)   
  
  
  point_alpha <- 0.2
  bgd_line_alpha <- point_alpha
  arrow_b_size <- 10
  max_point_size <- 7
  min_point_size <- 1
  bgd_line_width <- min_point_size + 1
  point_stroke = 0.25
  breaks <- seq(0, 20, 0.25)
  
  b_size = 12
  theme_set(theme_bw(base_size = b_size))
  
  
  # plot
  p <- smooth_data %>%   
    arrange(NormTime, Vowel) %>% 
    ggplot(aes(x = F2, y = F1, 
               color = Vowel, fill = Vowel)) +
    geom_path(linewidth = bgd_line_width, alpha = bgd_line_alpha) +
    geom_point(aes(size = veloc), color = "black", alpha = point_alpha, shape = 21, stroke = point_stroke) +
    geom_point(data = arrow_ends, aes(size = veloc - arrow_b_size), shape = 24,
               color = "black", alpha = 0.75, stroke = point_stroke) +
    scale_x_reverse("F2 (logHz)", breaks = breaks) +
    scale_y_reverse("F1 (logHz)", breaks = breaks) +
    scale_fill_manual(values = colors, 
                      guide = "none"
    ) +
    # guides(fill = guide_legend(override.aes = list(size = 5))) +
    scale_color_manual(values = colors, 
                       guide = "none"
    ) +
    scale_size_continuous(range = c(max_point_size, min_point_size), limits = c(-arrow_b_size, NA), guide = "none")
  
  # save
  ggsave(p,
         filename = "American_English_vowels_b_d- small.png",
         path = dir1,
         dpi = "retina", width = unit(6, "inches"), height = unit(6, "inches"))
  
  
## II  Vowel space heart animation w/ audio -----------------------------
  
  
  dir2 <- "Vowelspace heart animation with audio/"
    
  
  heart <- paste0(dir2, "SFlego_vowelspace_heart_formants_hand_corrected.csv") %>% 
    read.csv(header = TRUE)
  
  
  # what is this data? (can change variables to 'F1' and 'F2' to see what
  # they looked like before hand correction)
  heart %>% 
    ggplot(aes(x = IntervalTime)) +
    geom_line(aes(y = F1c)) +
    geom_line(aes(y = F2c)) +
    scale_x_continuous("Time (s)", breaks = seq(0, 15, 1)) +
    scale_y_continuous("F1 & F2 Frequency (Hz)", breaks = seq(0, 4000, 500))
  
  
  
  # start and end of the vocalization
  beginning <- 1.589086
  end <- 12.264709
  
  
  
  
  
                   
  # create a smooth contour over the vocalization with GAMs
  {
    
    pred_F1 <-  gam(F1c ~ s(IntervalTime, k = 100), data = heart)
    pred_F2 <-  gam(F2c ~ s(IntervalTime, k = 100), data = heart)               
    
    
    smooth_heart <- heart %>% 
      group_by(Filename) %>% 
      reframe(Duration = mean(unique(Duration)),
              IntervalTime = seq(0, Duration, 0.001)) %>% 
      mutate(NormTime = IntervalTime/Duration, .after = IntervalTime)
    
    
    smooth_heart$sF1 <- predict(pred_F1, newdata = smooth_heart, type = "response")
    smooth_heart$sF2 <- predict(pred_F2, newdata = smooth_heart, type = "response")
    
    
    smooth_heart <- smooth_heart %>%
      mutate(
        sF1 = ifelse(IntervalTime < beginning | IntervalTime > end, NA, round(sF1)),
        sF2 = ifelse(IntervalTime < beginning | IntervalTime > end, NA, round(sF2)))
    
    
    rm(pred_F1, pred_F2)
    
  }
  
  
  
  smooth_heart <- smooth_heart %>% 
    filter(IntervalTime > beginning & IntervalTime < end) %>% 
    group_by(Filename) %>% 
    mutate(F1_veloc = predict(smooth.spline(x = IntervalTime, y = sF1), deriv = 1)$y ) %>%
    mutate(F2_veloc = predict(smooth.spline(x = IntervalTime, y = sF2), deriv = 1)$y ) %>% 
    ungroup() %>% 
    mutate(VS_veloc = sqrt(F1_veloc^2 + F2_veloc^2)) %>% 
    right_join(smooth_heart) %>% 
    arrange(NormTime)
  
## II.A  plotting parameters that work well for 6" static image -----------------------------
  

  point_alpha <- 0.2
  bgd_line_alpha <- point_alpha
  max_point_size <- 7
  min_point_size <- 1
  bgd_line_width <- min_point_size + 1
  point_stroke <- 0.35
  
  b_size = 12
  theme_set(theme_bw(base_size = b_size))
  
  
  p <-  smooth_heart %>%   
    mutate(IntervalTime2 = IntervalTime) %>% 
    ggplot(aes(x = sF2, y = sF1)) +
    geom_path(linewidth = bgd_line_width, color = "#C51104", alpha = bgd_line_alpha) +
    geom_point(aes(size = VS_veloc, group = IntervalTime), fill = "#C51104",
               color = "#290000", alpha = point_alpha, shape = 21, stroke = point_stroke) +
    scale_x_reverse("F2 (Hz)") +
    scale_y_reverse("F1 (Hz)") +
    scale_size_continuous(range = c(max_point_size, min_point_size), limits = c(0, NA), guide = "none") +
    theme(legend.position = "none")
  
  
  ggsave(p, 
         filename = "SFlego_vowelspace_heart_static_image.png",
         path = dir2,
         dpi = "retina", width = unit(6, "inches"), height = unit(6, "inches"))
  
  
## II.B Create animated versions w/ and w/o synced audio -------------------
  
  
  
  anim <- animate(p + 
                  transition_reveal(IntervalTime),
                  # nframes = round(nrow(smooth)/10), # can use this instead of 'fps'
                  duration = max(smooth_heart$IntervalTime),
                  fps = 1.000 / 0.010, # every frame takes 10 ms (fps of 100)
                  renderer = av_renderer("SFlego_vowelspace_heart_animated_no_audio.mp4", dir = dir2),
                  height = 6, width = 6, units = "in", res = 320) 
  
  
  
  
  ## Add .wav audio to .mp4 animation 
  av_encode_video(
    input = paste0(dir2, "SFlego_vowelspace_heart_animated_no_audio.mp4"),   # the video input
    audio = paste0(dir2,"SFlego_vowelspace_heart_audio.wav"),               # the audio input
    output = paste0(dir2,"SFlego_vowelspace_heart_animated_w_audio.mp4"),   # the combined output
    framerate = 1.000 / 0.010 # same frame rate as above when first generated
  )
  
  
  
  
  
  
  
  
  
  
  
  
  

  