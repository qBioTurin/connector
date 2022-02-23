---
title: "CONNECTOR - Home"
layout: splash
header:
  overlay_filter: "0.7"
  overlay_image: /assets/images/HomeLogo.png
  actions:
    - label: "Learn More"
      url: "/framework/"
excerpt: "CONNECTOR"
intro: 
  - excerpt: '**CONNECOTR** represents ..........'
feature_row:
  - image_path: /assets/images/COVID/COVIDmodel.png
    alt: "placeholder image 2"
    title: "COVID-19"
    excerpt: "Investigation of the COVID-19 diffusion in the Piedmonnt region"
    url: "/covid19/"
    btn_label: "Read More"
    btn_class: "btn--primary"  
  - image_path: /assets/images/Pertussis/PertussisModel.png
    alt: "placeholder image 2"
    title: "Pertussis"
    excerpt: "Investigation of the pertussis epidemiology in Italy"
    url: "/Pertussis/"
    btn_label: "Read More"
    btn_class: "btn--primary"  
  - image_path: /assets/images/MS/MSmodel.jpg
    alt: "placeholder image 2"
    title: "Multiple Sclerosis"
    excerpt: "Analysis of the immune response in Multiple Sclerosis given specific treatments"
    url: "/ms/"
    btn_label: "Read More"
    btn_class: "btn--primary"  
feature_rowFramework:
  - image_path: /assets/images/Framework.png
    alt: "placeholder Framework"
    title: "Framework"
    excerpt: 'New general modeling framework for the analysis of epidemiological and biological systems, which exploits Petri Net graphical formalism, R environment, and Docker containerization to derive a tool easily accessible by any researcher even without advanced mathematical and computational skills.'
    url: /framework/
    btn_label: "Read More"
    btn_class: "btn--primary"
---

{% include feature_row id="intro" type="center" %}

{% include feature_row id="feature_rowFramework" type="left" %}

{% include feature_row %}