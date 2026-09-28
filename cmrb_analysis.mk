# ==============================================================================
# 累积代谢风险负担与肾功能快速下降：完整分析管线（修正版）
# 对应文章：Cumulative Metabolic Risk Burden and Rapid Decline in Renal Function:
#           A Longitudinal Cohort Study with Null Findings
#
# 相对 cmrb1.Rmd 的关键修正：
#   [FIX-1] 数据文件实为 25 列且首行是表头；原代码用 col_names=FALSE 并截取前 24 列，
#           导致全部变量整体左移一位：血小板被当作 eGFR、BMI 被当作收缩压、
#           收缩压被当作 BMI、总胆固醇被当作 HDL-C、既往史文本被当作空腹血糖(全变 NA)，
#           而真正的 EGFR 列被整列丢弃。现按表头正确读取。
#   [FIX-2] 同人同年重复体检（409 个人年）未去重，会使斜率回归重复计权。现按信息完整度去重。
#   [FIX-3] 原代码未实施文章 Methods 的入选标准（age>=18、基线 eGFR>=60、排除 ESRD、
#           关键变量缺失>2 次随访者排除）。现完整实施并输出流程图计数。
#   [FIX-4] 年度代谢得分原用 rowSums(na.rm=TRUE)，把缺失组分静默当作"正常"，
#           系统性低估负担。现缺失组分不参与计分，该年度得分记为 NA。
#   [FIX-5] 未清洗生理不可能的录入错误（如 eGFR=4412 对应肌酐 2.7 µmol/L）。现按区间置缺。
#   [ADD-1] 补齐文章要求但原代码缺失的：Model 1/2/3 汇总表、Wald χ²、四分位森林图、
#           NRI 与 IDI、MICE 多重插补、性别与年龄交互检验、eGFR 轨迹图、组分患病率图。
# ==============================================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(stringr)
  library(pROC); library(ggplot2); library(tableone); library(broom)
  library(mice); library(scales); library(patchwork)
})

setwd("D:/rProject/cmrb")
dir.create("figures", showWarnings = FALSE)
dir.create("tables",  showWarnings = FALSE)

LOG <- file("results_summary.txt", open = "wt", encoding = "UTF-8")
say <- function(...) { msg <- paste0(...); cat(msg, "\n"); cat(msg, "\n", file = LOG) }
say_sep <- function(t) { s <- paste0("\n", strrep("=", 74), "\n", t, "\n", strrep("=", 74));
                         cat(s, "\n"); cat(s, "\n", file = LOG) }

# 出版级图形主题
theme_pub <- function(base = 11) {
  theme_classic(base_size = base) +
    theme(plot.title = element_text(face = "bold", size = base + 1, hjust = 0),
          plot.subtitle = element_text(size = base - 1, colour = "grey30"),
          axis.title = element_text(face = "bold"),
          legend.position = "bottom",
          legend.title = element_text(face = "bold"),
          strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text = element_text(face = "bold"))
}
save_fig <- function(p, name, w, h) {
  ggsave(file.path("figures", name), p, width = w, height = h, dpi = 300, bg = "white")
  say("  [图] ", name, "  (", w, "x", h, " in)")
}

# ==============================================================================
# 1. 读取数据（按表头读取，映射为规范英文变量名）
# ==============================================================================
say_sep("1. 数据读取")

file_paths <- c("data_16.xlsx","data_17.xlsx","data_18.xlsx","data_19.xlsx","data_20.xlsx")
years      <- 2016:2020

# 真实表头 -> 规范变量名
col_map <- c(
  "SFZH" = "id", "XB" = "sex", "age" = "age_visit", "体重" = "weight",
  "体重指数" = "bmi", "本次收缩压(mmHg)" = "sbp", "本次舒张压(mmHg)" = "dbp",
  "饮酒史" = "alcohol", "吸烟史" = "smoking", "既往史" = "hx",
  "空腹血糖(GLU)" = "glu", "甘油三酯(TG)" = "tg", "总胆固醇(TC)" = "tc",
  "高密度脂蛋白(HDL_C)" = "hdl", "尿酸(UA)" = "ua", "尿蛋白定性" = "upro",
  "隐血(BLD)" = "bld", "低密度脂蛋白(LDL)" = "ldl", "肌酐(Cre)" = "cre",
  "总蛋白" = "tp", "白蛋白" = "alb", "白细胞(WBC)" = "wbc",
  "血红蛋白(HGB)" = "hgb", "血小板(PLT)" = "plt", "EGFR" = "egfr"
)

# 体检记录中用 "-" / 空白等占位符表示未测，需转为 NA
clean_num <- function(x) {
  s <- trimws(as.character(x))
  s[s %in% c("", "-", "--", "—", "－", "NA", "N/A", "null", "NULL", "无")] <- NA_character_
  suppressWarnings(as.numeric(s))
}

num_src <- c("age_visit","weight","bmi","sbp","dbp","glu","tg","tc","hdl","ua","ldl",
             "cre","tp","alb","wbc","hgb","plt","egfr")

raw_list <- list()
for (i in seq_along(file_paths)) {
  d <- read_excel(file_paths[i])
  stopifnot(all(names(col_map) %in% names(d)))
  d <- d[, names(col_map)]
  names(d) <- unname(col_map)
  for (v in num_src) d[[v]] <- clean_num(d[[v]])
  d$hx      <- trimws(as.character(d$hx))
  d$sex     <- trimws(as.character(d$sex))
  d$id      <- trimws(as.character(d$id))
  d$year    <- years[i]
  d$file    <- file_paths[i]
  raw_list[[i]] <- d
  say(sprintf("  %s : %6d 行 × %2d 列  (年份 %d)", file_paths[i], nrow(d), ncol(d) - 2, years[i]))
}
raw <- bind_rows(raw_list)
say(sprintf("\n  合并总记录数 = %d   唯一身份证号 = %d", nrow(raw), n_distinct(raw$id)))

# ==============================================================================
# 2. 生理不可能值置缺（录入错误清洗）  [FIX-5]
# ==============================================================================
say_sep("2. 异常值清洗（生理不可能值置为缺失）")

# 合理生理区间；区间外视为录入错误
ranges <- list(egfr=c(5,250), cre=c(10,1500), sbp=c(70,260), dbp=c(40,160),
               bmi=c(12,60), weight=c(25,250), glu=c(1.5,40), tg=c(0.05,30),
               tc=c(1,25), hdl=c(0.1,5), ldl=c(0.1,15), ua=c(30,2000),
               age_visit=c(0,110), plt=c(10,1000), hgb=c(30,250), wbc=c(0.5,100))
n_flagged <- 0
for (v in names(ranges)) {
  lo <- ranges[[v]][1]; hi <- ranges[[v]][2]
  bad <- !is.na(raw[[v]]) & (raw[[v]] < lo | raw[[v]] > hi)
  if (any(bad)) {
    say(sprintf("  %-10s 区间外 %5d 条  (保留区间 %g ~ %g)", v, sum(bad), lo, hi))
    raw[[v]][bad] <- NA_real_; n_flagged <- n_flagged + sum(bad)
  }
}
say(sprintf("  合计置缺 %d 个单元格", n_flagged))

# ==============================================================================
# 3. 同人同年重复体检去重  [FIX-2]
# ==============================================================================
say_sep("3. 同人同年重复记录去重")

key_vars <- c("egfr","sbp","dbp","glu","tg","hdl","bmi","cre")
raw$n_key_ok <- rowSums(!is.na(raw[key_vars]))
dup_pairs <- raw %>% count(id, year) %>% filter(n > 1)
say(sprintf("  重复的(人-年)组合 = %d，涉及 %d 人", nrow(dup_pairs), n_distinct(dup_pairs$id)))
# 保留关键变量最完整的一条；并列时保留文件中靠前者
dat <- raw %>%
  group_by(id, year) %>%
  slice(which.max(n_key_ok)[1]) %>%
  ungroup() %>%
  select(-n_key_ok)
say(sprintf("  去重后记录数 = %d，唯一人数 = %d", nrow(dat), n_distinct(dat$id)))
n_visits_all <- dat %>% group_by(id) %>% summarise(n_visits = n(), .groups = "drop")
say("  每人随访次数分布：")
print(table(n_visits_all$n_visits))
say(paste("   ", paste(names(table(n_visits_all$n_visits)), table(n_visits_all$n_visits),
                      sep="次=", collapse="  ")))

# ==============================================================================
# 4. 年龄与既往史衍生变量
# ==============================================================================
say_sep("4. 衍生变量（年龄、病史标记）")

# 基线年龄以 2016 年为参照，从身份证号提取出生日期（文章 Methods）
id2birth <- function(x) {
  x <- as.character(x)
  b <- ifelse(nchar(x) == 18, substr(x, 7, 14),
       ifelse(nchar(x) == 15, paste0("19", substr(x, 7, 12)), NA_character_))
  as.Date(b, format = "%Y%m%d")
}
birth <- id2birth(dat$id)
dat$age_base_2016 <- 2016 - as.numeric(format(birth, "%Y"))
# 身份证不可解析时回退到该年体检记录的年龄列换算至 2016
fallback <- is.na(dat$age_base_2016) & !is.na(dat$age_visit)
dat$age_base_2016[fallback] <- dat$age_visit[fallback] - (dat$year[fallback] - 2016)
say(sprintf("  年龄由身份证解析成功 %d 条，回退用体检年龄列 %d 条，仍缺失 %d 条",
            sum(!fallback & !is.na(dat$age_base_2016)), sum(fallback),
            sum(is.na(dat$age_base_2016))))
say(sprintf("  基线年龄(2016参照) 全体 mean=%.2f sd=%.2f range=%g~%g",
            mean(dat$age_base_2016, na.rm=TRUE), sd(dat$age_base_2016, na.rm=TRUE),
            min(dat$age_base_2016, na.rm=TRUE), max(dat$age_base_2016, na.rm=TRUE)))

# 病史标记：既往史为自由文本，按关键词识别；"无高血压服药治疗"等否定表述需排除
has_word <- function(txt, pat) {
  txt <- ifelse(is.na(txt), "", txt)
  # 去掉"无X"形式的否定短语，避免"无高血压服药治疗"被误判为高血压
  neg <- c("无高血压", "无糖尿", "无血糖", "无血脂", "无血压高", "血压偏低", "血压正常",
           "无高血脂", "不高", "正常")
  hit <- grepl(pat, txt)
  for (g in neg) hit <- hit & !grepl(g, txt)
  hit
}
dat$hx_htn <- as.integer(has_word(dat$hx, "高血压|血压高|血压偏高|血压不稳|血压临界|血压偶尔升高"))
dat$hx_dm  <- as.integer(has_word(dat$hx, "糖尿病|血糖高|高血糖|血糖偏高|血糖略高|血糖异常"))
dat$hx_esrd <- as.integer(grepl("尿毒症|透析|肾脏替代", ifelse(is.na(dat$hx), "", dat$hx)))

say(sprintf("  自报高血压史记录数 = %d (%.1f%%)", sum(dat$hx_htn), 100*mean(dat$hx_htn)))
say(sprintf("  自报糖尿病/高血糖史记录数 = %d (%.1f%%)", sum(dat$hx_dm), 100*mean(dat$hx_dm)))
say(sprintf("  自报 ESRD/透析记录数 = %d", sum(dat$hx_esrd)))

# ==============================================================================
# 5. 五项代谢异常定义（文章 Methods，依中国指南）
# ==============================================================================
say_sep("5. 代谢异常组分定义")

dat$abn_htn <- as.integer(dat$sbp >= 140 | dat$dbp >= 90 | dat$hx_htn == 1)
dat$abn_glu <- as.integer(dat$glu >= 6.1 | dat$hx_dm == 1)
dat$abn_tg  <- as.integer(dat$tg  >= 1.7)
dat$abn_hdl <- as.integer((dat$sex == "男" & dat$hdl < 1.0) |
                          (dat$sex == "女" & dat$hdl < 1.3))
dat$abn_ob  <- as.integer(dat$bmi >= 28)
abn_cols <- c("abn_htn","abn_glu","abn_tg","abn_hdl","abn_ob")

# [FIX-4] 缺失组分不计为 0：该组分记 NA，年度得分需五项齐全才有效
for (v in abn_cols) dat[[v]][is.na(dat[[v]])] <- NA_integer_
dat$n_abn_known <- rowSums(!is.na(dat[abn_cols]))
dat$met_score   <- ifelse(dat$n_abn_known == 5,
                          rowSums(dat[abn_cols], na.rm = TRUE), NA_integer_)

say("  各组分患病率（全部记录，分母为非缺失）：")
lbl <- c(abn_htn="高血压", abn_glu="高血糖", abn_tg="高甘油三酯",
         abn_hdl="低HDL-C", abn_ob="肥胖(BMI>=28)")
prev_tab <- data.frame(
  component = unname(lbl[abn_cols]), var = abn_cols,
  n_valid = sapply(abn_cols, function(v) sum(!is.na(dat[[v]]))),
  n_abn   = sapply(abn_cols, function(v) sum(dat[[v]] == 1, na.rm = TRUE)),
  stringsAsFactors = FALSE)
prev_tab$prev_pct <- 100 * prev_tab$n_abn / prev_tab$n_valid
print(prev_tab, row.names = FALSE)
say(sprintf("\n  年度代谢得分(0~5) 非缺失 %d 条，mean=%.2f sd=%.2f",
            sum(!is.na(dat$met_score)), mean(dat$met_score, na.rm=TRUE), sd(dat$met_score, na.rm=TRUE)))
say("  年度得分分布：")
say(paste("   ", paste(names(table(dat$met_score, useNA="ifany")),
                       table(dat$met_score, useNA="ifany"), sep="分=", collapse="  ")))

# ==============================================================================
# 6. 入选 / 排除流程（文章 Methods）  [FIX-3]
# ==============================================================================
say_sep("6. 研究人群构建（入选/排除流程）")

flow <- list()
flow[["源人群：2016-2020 年度体检记录"]] <- n_distinct(dat$id)

# 标准1：基线年龄 >= 18 岁
keep_age <- dat %>% filter(!is.na(age_base_2016), age_base_2016 >= 18) %>% pull(id) %>% unique()
s1 <- setdiff(dat$id, keep_age)
flow[[sprintf("排除①基线年龄<18岁或年龄缺失 (n=%d)", length(s1))]] <- length(keep_age)
dat <- dat %>% filter(id %in% keep_age)

# 标准2：排除基线 ESRD / 肾脏替代治疗
esrd_ids <- dat %>% filter(hx_esrd == 1) %>% pull(id) %>% unique()
keep_esrd <- setdiff(unique(dat$id), esrd_ids)
flow[[sprintf("排除②基线ESRD或肾脏替代治疗 (n=%d)", length(esrd_ids))]] <- length(keep_esrd)
dat <- dat %>% filter(id %in% keep_esrd)

# 标准3：关键变量缺失的随访次数 > 2 者排除
crit <- c("egfr","sbp","dbp","glu","tg","hdl","bmi")
miss_cnt <- dat %>% group_by(id) %>%
  summarise(n_miss_visit = sum(rowSums(is.na(pick(all_of(crit)))) > 0),
            n_visits = n(), .groups = "drop")
excl_miss <- miss_cnt %>% filter(n_miss_visit > 2) %>% pull(id)
keep_miss <- setdiff(unique(dat$id), excl_miss)
flow[[sprintf("排除③关键变量缺失随访次数>2 (n=%d)", length(excl_miss))]] <- length(keep_miss)
dat <- dat %>% filter(id %in% keep_miss)
miss_cnt <- miss_cnt %>% filter(id %in% keep_miss)

# 标准4：至少完成 3 次年检
keep_n3 <- miss_cnt %>% filter(n_visits >= 3) %>% pull(id)
excl_n3 <- setdiff(unique(dat$id), keep_n3)
flow[[sprintf("排除④随访次数<3次 (n=%d)", length(excl_n3))]] <- length(keep_n3)
dat <- dat %>% filter(id %in% keep_n3)

# 标准5：基线（首次随访）eGFR >= 60
first_visit <- dat %>% group_by(id) %>% summarise(base_year = min(year), .groups = "drop")
base_egfr_tab <- dat %>% inner_join(first_visit, by = c("id","year" = "base_year")) %>%
  select(id, base_year = year, egfr_base = egfr, sbp_base = sbp, bmi_base = bmi,
         sex_base = sex, age_base = age_base_2016, met_base = met_score)
keep_egfr <- base_egfr_tab %>% filter(!is.na(egfr_base), egfr_base >= 60) %>% pull(id)
excl_egfr <- setdiff(unique(dat$id), keep_egfr)
flow[[sprintf("排除⑤基线eGFR<60或缺失 (n=%d)", length(excl_egfr))]] <- length(keep_egfr)
dat <- dat %>% filter(id %in% keep_egfr)
base_egfr_tab <- base_egfr_tab %>% filter(id %in% keep_egfr)

# 标准6：eGFR 斜率可估（>=2 个非缺失 eGFR 且年份有变异）
slope_tab <- dat %>% filter(!is.na(egfr)) %>% group_by(id) %>%
  summarise(n_egfr = n(),
            slope = if (n() >= 2 && sd(year) > 0) unname(coef(lm(egfr ~ year))[2]) else NA_real_,
            slope_se = if (n() >= 3 && sd(year) > 0)
              summary(lm(egfr ~ year))$coefficients[2, 2] else NA_real_,
            follow_years = max(year) - min(year) + 1,
            last_year = max(year), .groups = "drop")
keep_slope <- slope_tab %>% filter(!is.na(slope)) %>% pull(id)
excl_slope <- setdiff(unique(dat$id), keep_slope)
flow[[sprintf("排除⑥eGFR斜率不可估 (n=%d)", length(excl_slope))]] <- length(keep_slope)
dat <- dat %>% filter(id %in% keep_slope)
slope_tab <- slope_tab %>% filter(id %in% keep_slope)

for (nm in names(flow)) say(sprintf("  %-52s 剩余 %6d 人", nm, flow[[nm]]))
n_final_ids <- unique(dat$id)
say(sprintf("\n  >>> 分析人群（结局可评估）= %d 人，人年 = %d", length(n_final_ids), nrow(dat)))

# ==============================================================================
# 7. 暴露（累积负担）与结局（快速下降）
# ==============================================================================
say_sep("7. 暴露与结局变量")

burden <- dat %>% filter(!is.na(met_score)) %>% group_by(id) %>%
  summarise(total_burden = sum(met_score), n_scored_visits = n(),
            mean_burden  = mean(met_score), .groups = "drop")

ana <- base_egfr_tab %>%
  left_join(slope_tab, by = "id") %>%
  left_join(burden,    by = "id") %>%
  left_join(n_visits_all, by = "id") %>%
  rename(age = age_base, sex = sex_base, base_egfr = egfr_base,
         base_sbp = sbp_base, base_bmi = bmi_base, base_met = met_base) %>%
  mutate(sex = factor(sex, levels = c("男","女"), labels = c("Male","Female")),
         rapid_decline = as.integer(slope < -5),
         rapid_decline_alt = as.integer(slope < -3))

say(sprintf("  有效得分年度数：mean=%.2f  分布：%s",
            mean(ana$n_scored_visits, na.rm=TRUE),
            paste(names(table(ana$n_scored_visits)), table(ana$n_scored_visits),
                  sep="年=", collapse=" ")))
say(sprintf("  累积负担 total_burden：mean=%.2f sd=%.2f median=%.0f range=%g~%g",
            mean(ana$total_burden, na.rm=TRUE), sd(ana$total_burden, na.rm=TRUE),
            median(ana$total_burden, na.rm=TRUE),
            min(ana$total_burden, na.rm=TRUE), max(ana$total_burden, na.rm=TRUE)))
say(sprintf("  平均负担 mean_burden ：mean=%.2f sd=%.2f range=%.2f~%.2f",
            mean(ana$mean_burden, na.rm=TRUE), sd(ana$mean_burden, na.rm=TRUE),
            min(ana$mean_burden, na.rm=TRUE), max(ana$mean_burden, na.rm=TRUE)))
say(sprintf("  eGFR 斜率(mL/min/1.73m²/年)：mean=%.3f sd=%.3f median=%.3f range=%.2f~%.2f",
            mean(ana$slope, na.rm=TRUE), sd(ana$slope, na.rm=TRUE), median(ana$slope, na.rm=TRUE),
            min(ana$slope, na.rm=TRUE), max(ana$slope, na.rm=TRUE)))
say(sprintf("  随访年数：median=%d range=%d~%d；人年合计=%d",
            median(ana$follow_years), min(ana$follow_years), max(ana$follow_years),
            sum(ana$follow_years)))
say(sprintf("\n  快速下降(slope<-5)：%d 人 (%.1f%%)", sum(ana$rapid_decline), 100*mean(ana$rapid_decline)))
say(sprintf("  宽松阈值(slope<-3)：%d 人 (%.1f%%)",
            sum(ana$rapid_decline_alt), 100*mean(ana$rapid_decline_alt)))

# 建模数据集：主要分析用完整病例
model_vars <- c("rapid_decline","rapid_decline_alt","total_burden","mean_burden","age",
                "sex","base_egfr","base_sbp","base_bmi","base_met","slope","n_visits",
                "follow_years","n_scored_visits")
ana_full <- ana %>% select(id, all_of(model_vars))
dm <- ana_full %>% filter(if_all(c("rapid_decline","total_burden","mean_burden","age",
                                   "sex","base_egfr","base_sbp","base_bmi","base_met"),
                                 ~ !is.na(.)))
say(sprintf("\n  完整病例建模样本量 = %d（排除协变量缺失 %d 人）",
            nrow(dm), nrow(ana_full) - nrow(dm)))
say("  缺失情况（分析人群 ana_full）：")
miss_sum <- data.frame(var = model_vars,
                       n_missing = sapply(model_vars, function(v) sum(is.na(ana_full[[v]]))))
miss_sum$pct <- round(100 * miss_sum$n_missing / nrow(ana_full), 2)
print(miss_sum[miss_sum$n_missing > 0, ], row.names = FALSE)
if (all(miss_sum$n_missing == 0)) say("    （无缺失）")

# 标准化协变量（文章：Continuous covariates were standardised）
dm <- dm %>% mutate(across(c(total_burden, mean_burden, age, base_egfr, base_sbp, base_bmi),
                           ~ as.numeric(scale(.)), .names = "z_{col}"))

write.csv(dm, "analysis_data_corrected.csv", row.names = FALSE)
say("\n  [数据] 建模数据集已保存 analysis_data_corrected.csv")

# ==============================================================================
# 8. Figure 1：研究流程图
# ==============================================================================
say_sep("8. 图表生成")

flow_df <- data.frame(
  step = c("Source population:\ncommunity adults with annual\nexaminations 2016-2020",
           "Age >= 18 years at baseline",
           "No ESRD / renal replacement\ntherapy at enrolment",
           "Key variables missing on\n<= 2 follow-up visits",
           ">= 3 annual examinations",
           "Baseline eGFR >= 60\nmL/min/1.73 m2",
           "eGFR slope estimable\n(>= 2 valid eGFR values)",
           "Final analytic sample"),
  n = c(unname(flow[[1]]), unname(flow[[2]]), unname(flow[[3]]), unname(flow[[4]]),
        unname(flow[[5]]), unname(flow[[6]]), unname(flow[[7]]), nrow(ana_full)),
  stringsAsFactors = FALSE)
flow_df$excl_lab <- c("", paste0("Excluded n = ", unname(flow[[1]]) - unname(flow[[2]])),
                      paste0("Excluded n = ", unname(flow[[2]]) - unname(flow[[3]])),
                      paste0("Excluded n = ", unname(flow[[3]]) - unname(flow[[4]])),
                      paste0("Excluded n = ", unname(flow[[4]]) - unname(flow[[5]])),
                      paste0("Excluded n = ", unname(flow[[5]]) - unname(flow[[6]])),
                      paste0("Excluded n = ", unname(flow[[6]]) - unname(flow[[7]])), "")
flow_df$y <- rev(seq_len(nrow(flow_df)))
flow_df$is_final <- flow_df$step == "Final analytic sample"
flow_excl <- flow_df[flow_df$excl_lab != "", ]

p_flow <- ggplot(flow_df) +
  geom_rect(aes(xmin = 0.5, xmax = 3.5, ymin = y - 0.42, ymax = y + 0.42,
                fill = is_final), colour = "grey25", linewidth = 0.4) +
  geom_text(aes(x = 2, y = y, label = paste0(step, "\n(n = ", format(n, big.mark = ","), ")")),
            size = 3.1, lineheight = 0.92, colour = "black") +
  geom_segment(aes(x = 2, xend = 2, y = y + 0.58, yend = y + 0.42),
               data = flow_df[-1, ], arrow = arrow(length = unit(0.18, "cm")),
               colour = "grey35", linewidth = 0.5) +
  geom_segment(aes(x = 3.5, xend = 4.15, y = y, yend = y),
               data = flow_excl, colour = "grey45", linewidth = 0.4) +
  geom_text(aes(x = 4.25, y = y, label = excl_lab),
            data = flow_excl, hjust = 0, size = 2.9, colour = "grey20") +
  scale_fill_manual(values = c("FALSE" = "grey95", "TRUE" = "#cfe3f5"), guide = "none") +
  scale_x_continuous(limits = c(0.4, 6.6), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(title = "Figure 2. Study participant selection flow",
       subtitle = "Community-based longitudinal health examination registry, 2016-2020") +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0),
        plot.subtitle = element_text(colour = "grey30", hjust = 0),
        plot.margin = margin(8, 12, 8, 8))
save_fig(p_flow, "fig2_flow_diagram.png", 8.6, 8.2)

# ==============================================================================
# 9. Table 1：基线特征（按快速下降分层）
# ==============================================================================
say_sep("9. Table 1 基线特征")

t1_vars <- c("age","sex","base_egfr","base_sbp","base_bmi","total_burden","mean_burden","base_met")
t1 <- CreateTableOne(vars = t1_vars, strata = "rapid_decline", data = dm,
                     factorVars = c("sex"), addOverall = TRUE, test = TRUE)
t1_out <- print(t1, showAllLevels = TRUE, formatOptions = list(big.mark = ",", digits = 2),
                printToggle = FALSE, quote = FALSE, noSpaces = TRUE)
write.csv(t1_out, "tables/table1_baseline_characteristics.csv")
say("  [表] tables/table1_baseline_characteristics.csv")
print(t1_out)

# 连续变量组间 t 检验 + 分类变量卡方（文章 Methods）
say("\n  组间比较（独立样本 t 检验 / 卡方检验）：")
cmp_rows <- list()
for (v in c("age","base_egfr","base_sbp","base_bmi","total_burden","mean_burden","base_met")) {
  a <- dm[[v]][dm$rapid_decline == 0]; b <- dm[[v]][dm$rapid_decline == 1]
  tt <- t.test(b, a)
  cmp_rows[[length(cmp_rows)+1]] <- data.frame(
    Variable = v,
    Non_rapid_mean_sd = sprintf("%.2f ± %.2f", mean(a), sd(a)),
    Rapid_mean_sd     = sprintf("%.2f ± %.2f", mean(b), sd(b)),
    Diff = sprintf("%.2f", mean(b) - mean(a)),
    t = sprintf("%.2f", tt$statistic), P_value = format.pval(tt$p.value, eps = 0.001, digits = 3),
    Method = "t-test", stringsAsFactors = FALSE)
}
ct <- table(dm$sex, dm$rapid_decline); chi <- chisq.test(ct)
cmp_rows[[length(cmp_rows)+1]] <- data.frame(
  Variable = "sex (Female, %)",
  Non_rapid_mean_sd = sprintf("%d (%.1f)", ct[1,1], 100*ct[1,1]/sum(ct[,1])),
  Rapid_mean_sd = sprintf("%d (%.1f)", ct[1,2], 100*ct[1,2]/sum(ct[,2])),
  Diff = "-", t = sprintf("chi2=%.2f", chi$statistic),
  P_value = format.pval(chi$p.value, eps = 0.001, digits = 3),
  Method = "chi-square", stringsAsFactors = FALSE)
cmp_tab <- bind_rows(cmp_rows)
write.csv(cmp_tab, "tables/table1_groupwise_tests.csv", row.names = FALSE)
print(cmp_tab, row.names = FALSE)

# ==============================================================================
# 10. Logistic 回归：Model 1 / 2 / 3 + 四分位 + 趋势检验
# ==============================================================================
say_sep("10. Logistic 回归分析")

fit_or <- function(m, term_label_map = NULL) {
  co <- summary(m)$coefficients
  ci <- suppressMessages(confint(m))          # profile likelihood CI
  out <- data.frame(
    term = rownames(co),
    beta = co[, 1], se = co[, 2],
    wald_chi2 = (co[, 1] / co[, 2])^2,
    OR = exp(co[, 1]), lower = exp(ci[, 1]), upper = exp(ci[, 2]),
    p = co[, 4], stringsAsFactors = FALSE)
  if (!is.null(term_label_map)) out$label <- ifelse(out$term %in% names(term_label_map),
                                                    term_label_map[out$term], out$term)
  out
}
lbl_map <- c("(Intercept)"="Intercept", total_burden="Cumulative burden (per 1 point)",
             age="Age (per 1 year)", sexFemale="Sex (female vs. male)",
             base_egfr="Baseline eGFR (per 1 mL/min/1.73 m2)",
             base_sbp="Baseline SBP (per 1 mmHg)", base_bmi="Baseline BMI (per 1 kg/m2)",
             base_met="Baseline metabolic abnormalities (per 1 count)",
             mean_burden="Mean burden (per 1 point)")

m1 <- glm(rapid_decline ~ total_burden, data = dm, family = binomial)
m2 <- glm(rapid_decline ~ total_burden + age + sex + base_egfr + base_sbp,
          data = dm, family = binomial)
m3 <- glm(rapid_decline ~ total_burden + age + sex + base_egfr + base_sbp +
            base_bmi + base_met, data = dm, family = binomial)

r1 <- fit_or(m1, lbl_map); r2 <- fit_or(m2, lbl_map); r3 <- fit_or(m3, lbl_map)
fmt_or <- function(r, term) {
  x <- r[r$term == term, ]
  sprintf("%.3f (%.3f-%.3f)", x$OR, x$lower, x$upper)
}
fmt_p <- function(r, term) format.pval(r$p[r$term == term], eps = 0.001, digits = 3)

model_sum <- data.frame(
  Variable = r2$label[r2$term != "(Intercept)"],
  Model1_crude_OR   = sapply(r2$term[r2$term != "(Intercept)"], function(t)
    if (t %in% r1$term) fmt_or(r1, t) else "-"),
  Model1_P = sapply(r2$term[r2$term != "(Intercept)"], function(t)
    if (t %in% r1$term) fmt_p(r1, t) else "-"),
  Model2_adjusted_OR = sapply(r2$term[r2$term != "(Intercept)"], function(t) fmt_or(r2, t)),
  Model2_P = sapply(r2$term[r2$term != "(Intercept)"], function(t) fmt_p(r2, t)),
  Model3_fully_adjusted_OR = sapply(r2$term[r2$term != "(Intercept)"], function(t)
    if (t %in% r3$term) fmt_or(r3, t) else "-"),
  Model3_P = sapply(r2$term[r2$term != "(Intercept)"], function(t)
    if (t %in% r3$term) fmt_p(r3, t) else "-"),
  stringsAsFactors = FALSE, row.names = NULL)
write.csv(model_sum, "tables/table_models_1_2_3.csv", row.names = FALSE)
say("  [表] tables/table_models_1_2_3.csv")
print(model_sum, row.names = FALSE)

# Table 2：Model 2 详表（β、SE、Wald χ²、OR 95%CI、P）
tab2 <- r2 %>% filter(term != "(Intercept)") %>%
  transmute(Variable = label,
            `beta coefficient` = sprintf("%.4f", beta),
            SE = sprintf("%.4f", se),
            `Wald chi2` = sprintf("%.2f", wald_chi2),
            `OR (95% CI)` = sprintf("%.3f (%.3f-%.3f)", OR, lower, upper),
            `P value` = format.pval(p, eps = 0.001, digits = 3),
            check.names = FALSE)
write.csv(tab2, "tables/table2_logistic_model2.csv", row.names = FALSE)
say("\n  [表] tables/table2_logistic_model2.csv  (Model 2)")
print(tab2, row.names = FALSE)
say(sprintf("\n  Model 2 拟合优度：AIC=%.1f  伪R2(McFadden)=%.4f  n=%d  事件数=%d",
            AIC(m2), 1 - logLik(m2)/logLik(glm(rapid_decline ~ 1, data=dm, family=binomial)),
            nrow(dm), sum(dm$rapid_decline)))
# Hosmer-Lemeshow 检验（自行实现，避免额外依赖）
hl_test <- function(y, phat, g = 10) {
  brk <- quantile(phat, probs = seq(0, 1, length.out = g + 1), na.rm = TRUE)
  brk <- unique(brk); brk[1] <- -Inf; brk[length(brk)] <- Inf
  grp <- cut(phat, breaks = brk, include.lowest = TRUE)
  n_g  <- tapply(rep(1, length(y)), grp, sum)
  obs1 <- tapply(y, grp, sum); exp1 <- tapply(phat, grp, sum)
  stat <- sum((obs1 - exp1)^2 / (exp1 * (1 - exp1 / n_g)), na.rm = TRUE)
  dfree <- length(obs1) - 2
  list(statistic = stat, parameter = dfree,
       p.value = pchisq(stat, dfree, lower.tail = FALSE))
}
say("  Hosmer-Lemeshow 检验（Model 2 拟合优度）：")
hl <- hl_test(dm$rapid_decline, fitted(m2), g = 10)
say(sprintf("    chi2=%.3f, df=%d, P=%.4f  (P>0.05 表示拟合可接受)",
            hl$statistic, hl$parameter, hl$p.value))

# 标准化尺度（文章称连续协变量已标准化）
m2z <- glm(rapid_decline ~ z_total_burden + z_age + sex + z_base_egfr + z_base_sbp,
           data = dm, family = binomial)
rz <- fit_or(m2z)
say("\n  Model 2 标准化尺度（每增加 1 个 SD）：")
for (t in rz$term[rz$term != "(Intercept)"])
  say(sprintf("    %-20s OR=%.3f (%.3f-%.3f) p=%s", t,
              rz$OR[rz$term==t], rz$lower[rz$term==t], rz$upper[rz$term==t],
              format.pval(rz$p[rz$term==t], eps=0.001, digits=3)))

# 四分位数分组（离散暴露，用分位数断点并合并重复切点）
qs <- quantile(dm$total_burden, c(0, .25, .5, .75, 1), na.rm = TRUE)
say(sprintf("\n  累积负担四分位断点：Q1<=%.0f, Q2<=%.0f, Q3<=%.0f, Q4<=%.0f",
            qs[2], qs[3], qs[4], qs[5]))
uq <- unique(qs)
if (length(uq) == length(qs)) {
  dm$burden_q <- cut(dm$total_burden, breaks = qs, include.lowest = TRUE,
                     labels = c("Q1","Q2","Q3","Q4"))
} else {
  say("  注：分布高度离散导致分位数断点重复，改用等距/实际取值分组")
  dm$burden_q <- ntile(dm$total_burden, 4)
  dm$burden_q <- factor(dm$burden_q, labels = c("Q1","Q2","Q3","Q4"))
}
dm$burden_q <- relevel(factor(dm$burden_q), ref = "Q1")
q_tab <- dm %>% group_by(burden_q) %>%
  summarise(n = n(), events = sum(rapid_decline),
            range_burden = sprintf("%g-%g", min(total_burden), max(total_burden)),
            median_burden = median(total_burden),
            event_rate = sprintf("%.1f%%", 100*mean(rapid_decline)), .groups = "drop")
say("  各组样本量与事件率："); print(as.data.frame(q_tab), row.names = FALSE)
write.csv(q_tab, "tables/quartile_distribution.csv", row.names = FALSE)

mq <- glm(rapid_decline ~ burden_q + age + sex + base_egfr + base_sbp,
          data = dm, family = binomial)
rq <- fit_or(mq)
dm$q_num <- as.numeric(dm$burden_q)
mtrend <- glm(rapid_decline ~ q_num + age + sex + base_egfr + base_sbp,
              data = dm, family = binomial)
rt <- fit_or(mtrend)
p_trend <- rt$p[rt$term == "q_num"]
say(sprintf("  跨四分位趋势检验 P = %s", format.pval(p_trend, eps=0.001, digits=3)))

quart_tab <- data.frame(
  Quartile = c("Q1 (reference)", "Q2", "Q3", "Q4"),
  `Burden range` = q_tab$range_burden,
  n = q_tab$n, Events = q_tab$events, `Event rate` = q_tab$event_rate,
  `Adjusted OR (95% CI)` = c("1.000 (ref)",
    sapply(c("burden_qQ2","burden_qQ3","burden_qQ4"), function(t) fmt_or(rq, t))),
  `P value` = c("-", sapply(c("burden_qQ2","burden_qQ3","burden_qQ4"),
                            function(t) fmt_p(rq, t))),
  check.names = FALSE, stringsAsFactors = FALSE, row.names = NULL)
quart_tab <- rbind(quart_tab, data.frame(Quartile = "P for trend", `Burden range` = "-",
  n = NA, Events = NA, `Event rate` = "-",
  `Adjusted OR (95% CI)` = sprintf("%.3f (%.3f-%.3f) per quartile",
    rt$OR[rt$term=="q_num"], rt$lower[rt$term=="q_num"], rt$upper[rt$term=="q_num"]),
  `P value` = format.pval(p_trend, eps=0.001, digits=3), check.names = FALSE))
write.csv(quart_tab, "tables/table_quartile_dose_response.csv", row.names = FALSE)
say("\n  [表] tables/table_quartile_dose_response.csv")
print(quart_tab, row.names = FALSE)

# 未调整四分位（便于与调整版对照）
mq_crude <- glm(rapid_decline ~ burden_q, data = dm, family = binomial)
rq_crude <- fit_or(mq_crude)

# ==============================================================================
# 11. 敏感性分析：平均负担 / 阈值 -3 / MICE
# ==============================================================================
say_sep("11. 敏感性分析")

sens_rows <- list()
add_sens <- function(label, r, term, n, ev) {
  x <- r[r$term == term, ]
  sens_rows[[length(sens_rows)+1]] <<- data.frame(
    Strategy = label, `Sample size` = n, Events = ev,
    `OR (95% CI)` = sprintf("%.3f (%.3f-%.3f)", x$OR, x$lower, x$upper),
    `P value` = format.pval(x$p, eps = 0.001, digits = 3),
    check.names = FALSE, stringsAsFactors = FALSE)
}
add_sens("Primary analysis (total burden, slope < -5)", r2, "total_burden",
         nrow(dm), sum(dm$rapid_decline))

ms <- glm(rapid_decline ~ mean_burden + age + sex + base_egfr + base_sbp,
          data = dm, family = binomial)
add_sens("Mean burden instead of total burden", fit_or(ms), "mean_burden",
         nrow(dm), sum(dm$rapid_decline))

ma <- glm(rapid_decline_alt ~ total_burden + age + sex + base_egfr + base_sbp,
          data = dm, family = binomial)
add_sens("Lenient threshold (slope < -3)", fit_or(ma), "total_burden",
         nrow(dm), sum(dm$rapid_decline_alt))

# 额外敏感性：Model 3 全调整；排除随访<5年者；限定基线eGFR 60-120
m3s <- fit_or(m3)
add_sens("Fully adjusted (Model 3 covariates)", m3s, "total_burden",
         nrow(dm), sum(dm$rapid_decline))
d5 <- dm %>% filter(follow_years == 5)
m5 <- glm(rapid_decline ~ total_burden + age + sex + base_egfr + base_sbp,
          data = d5, family = binomial)
add_sens("Restricted to 5 complete follow-up years", fit_or(m5), "total_burden",
         nrow(d5), sum(d5$rapid_decline))
dN <- dm %>% filter(base_egfr >= 60, base_egfr <= 120)
mN <- glm(rapid_decline ~ total_burden + age + sex + base_egfr + base_sbp,
          data = dN, family = binomial)
add_sens("Restricted to baseline eGFR 60-120 (no hyperfiltration)", fit_or(mN),
         "total_burden", nrow(dN), sum(dN$rapid_decline))

# MICE 多重插补（文章 Methods 要求）
say("\n  MICE 多重插补中（m=5, maxit=10, method=pmm）...")
imp_src <- ana_full %>% select(rapid_decline, total_burden, mean_burden, age, sex,
                               base_egfr, base_sbp, base_bmi, base_met) %>%
  mutate(rapid_decline = factor(rapid_decline, levels = c(0,1), labels = c("No","Yes")),
         sex = factor(sex))
n_missing_any <- sum(!complete.cases(imp_src))
say(sprintf("  插补前存在缺失的记录数 = %d / %d", n_missing_any, nrow(imp_src)))
if (n_missing_any > 0) {
  set.seed(20260921)
  imp <- mice(imp_src, m = 5, method = "pmm", maxit = 10, printFlag = FALSE, seed = 20260921)
  fit_imp <- with(imp, glm(rapid_decline ~ total_burden + age + sex + base_egfr + base_sbp,
                           family = binomial))
  pooled <- pool(fit_imp)
  ps <- summary(pooled, conf.int = TRUE)
  x <- ps[ps$term == "total_burden", ]
  sens_rows[[length(sens_rows)+1]] <- data.frame(
    Strategy = "Multiple imputation (MICE, m=5)", `Sample size` = nrow(imp_src),
    Events = sum(imp_src$rapid_decline == "Yes"),
    `OR (95% CI)` = sprintf("%.3f (%.3f-%.3f)", exp(x$estimate),
                            exp(x$`2.5 %`), exp(x$`97.5 %`)),
    `P value` = format.pval(x$p.value, eps = 0.001, digits = 3),
    check.names = FALSE, stringsAsFactors = FALSE)
  say("  MICE 全部变量合并估计：")
  ps_out <- ps %>% mutate(OR = exp(estimate), L95 = exp(`2.5 %`), U95 = exp(`97.5 %`)) %>%
    select(term, OR, L95, U95, p.value)
  print(as.data.frame(ps_out), row.names = FALSE, digits = 4)
  write.csv(as.data.frame(ps_out), "tables/sensitivity_mice_pooled.csv", row.names = FALSE)
} else {
  say("  完整病例无缺失，MICE 结果与主分析一致，按主分析值登记。")
  sens_rows[[length(sens_rows)+1]] <- data.frame(
    Strategy = "Multiple imputation (MICE) - no missing data, identical to primary",
    `Sample size` = nrow(dm), Events = sum(dm$rapid_decline),
    `OR (95% CI)` = fmt_or(r2, "total_burden"),
    `P value` = fmt_p(r2, "total_burden"), check.names = FALSE, stringsAsFactors = FALSE)
}

tab3 <- bind_rows(sens_rows)
write.csv(tab3, "tables/table3_sensitivity.csv", row.names = FALSE)
say("\n  [表] tables/table3_sensitivity.csv")
print(tab3, row.names = FALSE)

# ==============================================================================
# 12. ROC / AUC / DeLong / NRI / IDI
# ==============================================================================
say_sep("12. 预测性能：ROC、AUC、DeLong、NRI、IDI")

roc_burden <- roc(dm$rapid_decline, dm$total_burden, quiet = TRUE, ci = TRUE, boot.n = 2000)
pred_m2 <- predict(m2, type = "response")
roc_m2    <- roc(dm$rapid_decline, pred_m2, quiet = TRUE, ci = TRUE, boot.n = 2000)
# 核心模型（不含累积负担）：age + sex + base_egfr + base_sbp
m_core <- glm(rapid_decline ~ age + sex + base_egfr + base_sbp, data = dm, family = binomial)
pred_core <- predict(m_core, type = "response")
roc_core <- roc(dm$rapid_decline, pred_core, quiet = TRUE, ci = TRUE, boot.n = 2000)
roc_m3 <- roc(dm$rapid_decline, predict(m3, type = "response"), quiet = TRUE, ci = TRUE, boot.n = 2000)

ci_str <- function(r) sprintf("%.3f (%.3f-%.3f)", auc(r), ci(r)[1], ci(r)[3])
say(sprintf("  AUC 累积负担单独        = %s", ci_str(roc_burden)))
say(sprintf("  AUC 核心模型(age+sex+eGFR+SBP) = %s", ci_str(roc_core)))
say(sprintf("  AUC Model 2(核心+负担)  = %s", ci_str(roc_m2)))
say(sprintf("  AUC Model 3(全调整)     = %s", ci_str(roc_m3)))

dt1 <- roc.test(roc_burden, roc_m2, method = "delong")
dt2 <- roc.test(roc_core, roc_m2, method = "delong")
# pROC 1.19 的 roc.test 返回 htest：$estimate 为两条曲线的 AUC（非差值），
# 差值需自行相减；$conf.int 为差值的 CI；$statistic 为 Z。
dl_diff <- function(t) {
  e <- unname(t$estimate)
  if (length(e) >= 2) e[1] - e[2] else NA_real_
}
dl_z    <- function(t) unname(t$statistic[1])
dl_ci   <- function(t) if (!is.null(t$conf.int) && length(t$conf.int) >= 2)
                           c(unname(t$conf.int[1]), unname(t$conf.int[2])) else c(NA_real_, NA_real_)
say(sprintf("\n  DeLong 检验 负担 vs Model2 ：ΔAUC=%.3f, Z=%.3f, P=%s",
            dl_diff(dt1), dl_z(dt1), format.pval(dt1$p.value, eps=0.001, digits=3)))
say(sprintf("  DeLong 检验 核心 vs Model2 ：ΔAUC=%.3f, Z=%.3f, P=%s",
            dl_diff(dt2), dl_z(dt2), format.pval(dt2$p.value, eps=0.001, digits=3)))
say(sprintf("    核心 vs Model2 ΔAUC 95%%CI = %.3f ~ %.3f",
            dl_ci(dt2)[1], dl_ci(dt2)[2]))

# 连续 NRI 与 IDI：比较核心模型 vs 加入累积负担后的模型
roc_tbl <- data.frame(
  Model = c("Cumulative burden alone", "Core model (age, sex, baseline eGFR, baseline SBP)",
            "Model 2 (core + cumulative burden)", "Model 3 (fully adjusted)"),
  AUC = c(as.numeric(auc(roc_burden)), as.numeric(auc(roc_core)),
          as.numeric(auc(roc_m2)), as.numeric(auc(roc_m3))),
  `95% CI` = c(sprintf("%.3f-%.3f", ci(roc_burden)[1], ci(roc_burden)[3]),
               sprintf("%.3f-%.3f", ci(roc_core)[1], ci(roc_core)[3]),
               sprintf("%.3f-%.3f", ci(roc_m2)[1], ci(roc_m2)[3]),
               sprintf("%.3f-%.3f", ci(roc_m3)[1], ci(roc_m3)[3])),
  check.names = FALSE, stringsAsFactors = FALSE)
delong_rows <- data.frame(
  Comparison = c("Burden alone vs. Model 2", "Core model vs. Model 2"),
  `Delta AUC` = c(dl_diff(dt1), dl_diff(dt2)),
  Z = c(dl_z(dt1), dl_z(dt2)),
  `P value (DeLong)` = c(format.pval(dt1$p.value, eps=0.001, digits=3),
                         format.pval(dt2$p.value, eps=0.001, digits=3)),
  check.names = FALSE, stringsAsFactors = FALSE)
write.csv(roc_tbl, "tables/roc_auc_summary.csv", row.names = FALSE)
write.csv(delong_rows, "tables/delong_tests.csv", row.names = FALSE)
say("\n  [表] tables/roc_auc_summary.csv, tables/delong_tests.csv")
print(roc_tbl, row.names = FALSE)

# NRI / IDI（连续型）
nri_idi_calc <- function(y, p_base, p_new, nboot = 1000, seed = 20260921) {
  set.seed(seed)
  calc <- function(idx) {
    yy <- y[idx]; pb <- p_base[idx]; pn <- p_new[idx]
    ev <- yy == 1; nv <- yy == 0
    # 事件组上移 / 非事件组下移 => 正贡献
    up_e <- mean(pn[ev]  > pb[ev])  - mean(pn[ev]  < pb[ev])
    dn_n <- mean(pn[nv]  < pb[nv])  - mean(pn[nv]  > pb[nv])
    nri <- up_e + dn_n
    # IDI = (ISnew - ISold)，基于 Brier 分解
    idi <- (mean(pn[ev]) - mean(pn[nv])) - (mean(pb[ev]) - mean(pb[nv]))
    c(nri = nri, nri_event = up_e, nri_nonevent = dn_n, idi = idi)
  }
  point <- calc(seq_along(y))
  n <- length(y)
  boots <- replicate(nboot, calc(sample.int(n, n, replace = TRUE)))
  se <- apply(boots, 1, sd, na.rm = TRUE)
  data.frame(metric = names(point), estimate = point, se = se,
             lower = point - 1.96*se, upper = point + 1.96*se,
             z = point/se, p = 2*pnorm(-abs(point/se)), row.names = NULL)
}
y_bin <- as.integer(dm$rapid_decline == 1)
ni <- nri_idi_calc(y_bin, pred_core, pred_m2, nboot = 500)
say("\n  连续 NRI / IDI（核心模型 -> 加入累积负担，bootstrap 500 次）：")
print(ni, row.names = FALSE, digits = 4)
ni_out <- ni %>% transmute(Metric = metric, Estimate = sprintf("%.4f", estimate),
  `95% CI` = sprintf("%.4f to %.4f", lower, upper),
  `P value` = format.pval(p, eps = 0.001, digits = 3), check.names = FALSE)
write.csv(ni_out, "tables/table_nri_idi.csv", row.names = FALSE)
say("  [表] tables/table_nri_idi.csv")

# ==============================================================================
# 13. Figure 2：ROC 曲线（文章 Figure 1）
# ==============================================================================
roc_df <- bind_rows(
  data.frame(fpr = 1 - roc_burden$specificities, tpr = roc_burden$sensitivities,
             Model = sprintf("Cumulative burden alone\nAUC = %.3f (%.3f-%.3f)",
                             auc(roc_burden), ci(roc_burden)[1], ci(roc_burden)[3])),
  data.frame(fpr = 1 - roc_m2$specificities, tpr = roc_m2$sensitivities,
             Model = sprintf("Multivariable model\n(age, sex, baseline eGFR, baseline SBP)\nAUC = %.3f (%.3f-%.3f)",
                             auc(roc_m2), ci(roc_m2)[1], ci(roc_m2)[3])),
  data.frame(fpr = 1 - roc_core$specificities, tpr = roc_core$sensitivities,
             Model = sprintf("Core model without burden\n(age, sex, baseline eGFR, baseline SBP)\nAUC = %.3f (%.3f-%.3f)",
                             auc(roc_core), ci(roc_core)[1], ci(roc_core)[3])))
roc_df$Model <- factor(roc_df$Model, levels = unique(roc_df$Model))

p_roc <- ggplot(roc_df, aes(x = fpr, y = tpr, colour = Model, linetype = Model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey50") +
  geom_path(linewidth = 0.9) +
  scale_colour_manual(values = c("#d1495b", "#2c6fbb", "#5aa469")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash")) +
  coord_equal(xlim = c(0,1), ylim = c(0,1)) +
  labs(x = "1 - Specificity (false positive rate)",
       y = "Sensitivity (true positive rate)",
       colour = NULL, linetype = NULL,
       title = "Figure 1. ROC curves for prediction of rapid kidney function decline",
       subtitle = sprintf("n = %s (%.1f%% events). DeLong P, core vs. burden-added = %s",
                          format(nrow(dm), big.mark=","), 100*mean(dm$rapid_decline),
                          format.pval(dt2$p.value, eps=0.001, digits=3))) +
  theme_pub(11) +
  theme(legend.text = element_text(size = 8.6, lineheight = 0.9),
        plot.subtitle = element_text(size = 8.6))
save_fig(p_roc, "fig1_roc_curves.png", 8.4, 6.4)

# ==============================================================================
# 14. Figure 3：累积负担分布 + 剂量-反应关系
# ==============================================================================
dd <- dm %>% group_by(total_burden) %>%
  summarise(n = n(), events = sum(rapid_decline), rate = mean(rapid_decline), .groups="drop")
n_max_dd <- max(dd$n)
p_dist <- ggplot(dd, aes(x = total_burden)) +
  geom_col(aes(y = n), fill = "#8fb8de", colour = "#4a7ba7", linewidth = 0.2) +
  geom_line(aes(y = rate * n_max_dd), colour = "#d1495b", linewidth = 0.9) +
  geom_point(aes(y = rate * n_max_dd), colour = "#d1495b", size = 1.8) +
  scale_y_continuous(labels = label_comma(),
                     sec.axis = sec_axis(~ . / n_max_dd,
                                         labels = percent_format(accuracy = 1),
                                         name = "Rapid decline rate")) +
  labs(x = "Cumulative metabolic burden score (5-year total, range 0-25)",
       y = "Number of participants",
       title = "Figure 3A. Distribution of cumulative burden",
       subtitle = sprintf("mean = %.2f, SD = %.2f, median = %.0f, range = %g-%g",
                          mean(dm$total_burden), sd(dm$total_burden),
                          median(dm$total_burden),
                          min(dm$total_burden), max(dm$total_burden))) +
  theme_pub(10)

p_dose <- ggplot(dm, aes(x = total_burden, y = rapid_decline)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"),
              se = TRUE, colour = "#2c6fbb", fill = "#bcd3ec", linewidth = 1) +
  geom_point(data = dd, aes(x = total_burden, y = rate), size = 1.6, colour = "grey25") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Cumulative metabolic burden score (5-year total)",
       y = "Probability of rapid decline",
       title = "Figure 3B. Dose-response relationship (logistic, with 95% CI)",
       subtitle = sprintf("Dots = observed event rate per score level; P for trend across quartiles = %s",
                          format.pval(p_trend, eps = 0.001, digits = 3))) +
  theme_pub(10)

save_fig(p_dist / p_dose, "fig3_burden_distribution_dose_response.png", 7.6, 9.0)

# ==============================================================================
# 15. Figure 4：四分位 OR 森林图
# ==============================================================================
forest_df <- data.frame(
  group = c("Q1 (reference)", "Q2", "Q3", "Q4"),
  or = c(1, rq$OR[rq$term=="burden_qQ2"], rq$OR[rq$term=="burden_qQ3"], rq$OR[rq$term=="burden_qQ4"]),
  lo = c(1, rq$lower[rq$term=="burden_qQ2"], rq$lower[rq$term=="burden_qQ3"], rq$lower[rq$term=="burden_qQ4"]),
  hi = c(1, rq$upper[rq$term=="burden_qQ2"], rq$upper[rq$term=="burden_qQ3"], rq$upper[rq$term=="burden_qQ4"]),
  # Q1 为参照组，粗 OR 按定义固定为 1（不可用截距项代替）
  or_c = c(1, rq_crude$OR[rq_crude$term=="burden_qQ2"],
              rq_crude$OR[rq_crude$term=="burden_qQ3"],
              rq_crude$OR[rq_crude$term=="burden_qQ4"]),
  lo_c = c(1, rq_crude$lower[rq_crude$term=="burden_qQ2"],
              rq_crude$lower[rq_crude$term=="burden_qQ3"],
              rq_crude$lower[rq_crude$term=="burden_qQ4"]),
  hi_c = c(1, rq_crude$upper[rq_crude$term=="burden_qQ2"],
              rq_crude$upper[rq_crude$term=="burden_qQ3"],
              rq_crude$upper[rq_crude$term=="burden_qQ4"]),
  n = q_tab$n, ev = q_tab$events, stringsAsFactors = FALSE)
forest_df$group <- factor(forest_df$group, levels = rev(forest_df$group))
forest_long <- bind_rows(
  forest_df %>% transmute(group, est = or,   lo,    hi,    n, Model = "Adjusted (Model 2)"),
  forest_df %>% transmute(group, est = or_c, lo = lo_c, hi = hi_c, n, Model = "Crude"))
forest_long$Model <- factor(forest_long$Model, levels = c("Crude","Adjusted (Model 2)"))

p_forest_q <- ggplot(forest_long, aes(x = est, y = group, colour = Model)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi), height = 0.18, linewidth = 0.7,
                orientation = "y", position = position_dodge(width = 0.42)) +
  geom_point(aes(size = n), position = position_dodge(width = 0.42)) +
  scale_colour_manual(values = c("Crude" = "grey45", "Adjusted (Model 2)" = "#2c6fbb")) +
  scale_size_continuous(range = c(2.2, 5.0), name = "n") +
  scale_x_continuous(limits = c(min(forest_long$lo)*0.93, max(forest_long$hi)*1.07)) +
  labs(x = "Odds ratio for rapid kidney function decline (95% CI)", y = NULL,
       colour = NULL,
       title = "Figure 4. Rapid decline risk by quartile of cumulative metabolic burden",
       subtitle = sprintf("Q1 = reference. Adjusted for age, sex, baseline eGFR and baseline SBP. P for trend = %s",
                          format.pval(p_trend, eps = 0.001, digits = 3))) +
  theme_pub(10.5)
save_fig(p_forest_q, "fig4_forest_quartiles.png", 8.0, 4.6)

# ==============================================================================
# 16. Figure 5：敏感性分析森林图
# ==============================================================================
parse_or <- function(s) {
  m <- str_extract(s, "[0-9.]+ \\([0-9.]+-[0-9.]+\\)")
  num <- as.numeric(str_extract_all(m, "[0-9.]+")[[1]])
  num
}
sens_plot <- tab3 %>%
  mutate(est = sapply(`OR (95% CI)`, function(s) parse_or(s)[1]),
         lo  = sapply(`OR (95% CI)`, function(s) parse_or(s)[2]),
         hi  = sapply(`OR (95% CI)`, function(s) parse_or(s)[3]),
         lab = str_wrap(Strategy, 52))
sens_plot$lab <- factor(sens_plot$lab, levels = rev(sens_plot$lab))
p_sens <- ggplot(sens_plot, aes(x = est, y = lab)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi), height = 0.2, colour = "#4a7ba7",
                linewidth = 0.7, orientation = "y") +
  geom_point(colour = "#2c6fbb", size = 2.8) +
  geom_text(aes(x = hi, label = sprintf("%.3f (%.3f-%.3f)", est, lo, hi)),
            hjust = -0.12, size = 2.9, colour = "grey20") +
  scale_x_continuous(limits = c(min(sens_plot$lo)*0.97, max(sens_plot$hi)*1.30)) +
  labs(x = "Odds ratio per 1-point increment in cumulative burden (95% CI)", y = NULL,
       title = "Figure 5. Sensitivity analyses",
       subtitle = "All models adjusted for age, sex, baseline eGFR and baseline SBP unless stated otherwise") +
  theme_pub(10) + theme(axis.text.y = element_text(size = 8.4, lineheight = 0.9))
save_fig(p_sens, "fig5_forest_sensitivity.png", 9.6, 5.0)

# ==============================================================================
# 17. 亚组分析：性别分层 + 交互；年龄交互  Table 4 / Figure 6
# ==============================================================================
say_sep("13. 亚组分析与交互检验")

fit_sub <- function(d) {
  m <- glm(rapid_decline ~ total_burden + age + base_egfr + base_sbp,
           data = d, family = binomial)
  fit_or(m)
}
sub_rows <- list()
for (g in c("Male","Female")) {
  dg <- dm %>% filter(sex == g)
  if (nrow(dg) < 50) next
  rr <- fit_sub(dg); x <- rr[rr$term == "total_burden", ]
  sub_rows[[length(sub_rows)+1]] <- data.frame(
    Subgroup = g, n = nrow(dg), Events = sum(dg$rapid_decline),
    OR = x$OR, lower = x$lower, upper = x$upper, p = x$p, stringsAsFactors = FALSE)
  say(sprintf("  %-8s n=%5d 事件=%4d  OR=%.3f (%.3f-%.3f) P=%s", g, nrow(dg),
              sum(dg$rapid_decline), x$OR, x$lower, x$upper,
              format.pval(x$p, eps=0.001, digits=3)))
}
# 性别交互
m_int_sex <- glm(rapid_decline ~ total_burden * sex + age + base_egfr + base_sbp,
                 data = dm, family = binomial)
co_int <- summary(m_int_sex)$coefficients
int_term <- grep("^total_burden:sex", rownames(co_int), value = TRUE)
p_int_sex <- co_int[int_term[1], 4]
say(sprintf("  性别交互 P_interaction = %s", format.pval(p_int_sex, eps=0.001, digits=3)))

# 年龄交互 + 年龄分层
m_int_age <- glm(rapid_decline ~ total_burden * age + sex + base_egfr + base_sbp,
                 data = dm, family = binomial)
co_age <- summary(m_int_age)$coefficients
p_int_age <- co_age["total_burden:age", 4]
say(sprintf("  年龄交互 P_interaction = %s", format.pval(p_int_age, eps=0.001, digits=3)))

dm$age_grp <- cut(dm$age, breaks = c(-Inf, 45, 60, Inf),
                  labels = c("<45 years", "45-60 years", ">60 years"))
for (g in levels(dm$age_grp)) {
  dg <- dm %>% filter(age_grp == g)
  if (nrow(dg) < 50) next
  rr <- fit_sub(dg); x <- rr[rr$term == "total_burden", ]
  sub_rows[[length(sub_rows)+1]] <- data.frame(
    Subgroup = g, n = nrow(dg), Events = sum(dg$rapid_decline),
    OR = x$OR, lower = x$lower, upper = x$upper, p = x$p, stringsAsFactors = FALSE)
  say(sprintf("  %-12s n=%5d 事件=%4d  OR=%.3f (%.3f-%.3f) P=%s", g, nrow(dg),
              sum(dg$rapid_decline), x$OR, x$lower, x$upper,
              format.pval(x$p, eps=0.001, digits=3)))
}
# 基线 eGFR 是否高滤过（文章讨论重点）
dm$hf_grp <- ifelse(dm$base_egfr > 120, "Baseline eGFR > 120\n(hyperfiltration)",
                    "Baseline eGFR 60-120")
for (g in unique(dm$hf_grp)) {
  dg <- dm %>% filter(hf_grp == g)
  if (nrow(dg) < 50) next
  rr <- fit_sub(dg); x <- rr[rr$term == "total_burden", ]
  sub_rows[[length(sub_rows)+1]] <- data.frame(
    Subgroup = g, n = nrow(dg), Events = sum(dg$rapid_decline),
    OR = x$OR, lower = x$lower, upper = x$upper, p = x$p, stringsAsFactors = FALSE)
  say(sprintf("  %-28s n=%5d 事件=%4d  OR=%.3f (%.3f-%.3f) P=%s",
              gsub("\n"," ",g), nrow(dg), sum(dg$rapid_decline), x$OR, x$lower, x$upper,
              format.pval(x$p, eps=0.001, digits=3)))
}
sub_tab <- bind_rows(sub_rows) %>%
  mutate(`P value` = format.pval(p, eps = 0.001, digits = 3),
         `OR (95% CI)` = sprintf("%.3f (%.3f-%.3f)", OR, lower, upper))
sub_tab$P_interaction <- c(format.pval(p_int_sex, eps=0.001, digits=3), "",
                           rep("", nrow(sub_tab) - 2))
tab4 <- sub_tab %>% select(Subgroup, n, Events, `OR (95% CI)`, `P value`, P_interaction)
write.csv(tab4, "tables/table4_subgroup.csv", row.names = FALSE)
say("\n  [表] tables/table4_subgroup.csv"); print(tab4, row.names = FALSE)

sub_plot <- sub_tab %>% mutate(lab = factor(str_wrap(gsub("\n"," ",Subgroup), 34),
                                            levels = rev(str_wrap(gsub("\n"," ",Subgroup), 34))))
p_sub <- ggplot(sub_plot, aes(x = OR, y = lab)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = lower, xmax = upper), height = 0.2, colour = "#5aa469",
                linewidth = 0.7, orientation = "y") +
  geom_point(colour = "#2f7d43", size = 2.8) +
  geom_text(aes(x = upper, label = sprintf("%.3f (%.3f-%.3f)", OR, lower, upper)),
            hjust = -0.12, size = 2.9, colour = "grey20") +
  scale_x_continuous(limits = c(min(sub_plot$lower)*0.96, max(sub_plot$upper)*1.34)) +
  labs(x = "Odds ratio per 1-point increment in cumulative burden (95% CI)", y = NULL,
       title = "Figure 6. Subgroup analyses",
       subtitle = sprintf("Adjusted for age, baseline eGFR and baseline SBP within strata. P_interaction: sex = %s, age = %s",
                          format.pval(p_int_sex, eps=0.001, digits=3),
                          format.pval(p_int_age, eps=0.001, digits=3))) +
  theme_pub(10) + theme(axis.text.y = element_text(size = 8.6, lineheight = 0.9))
save_fig(p_sub, "fig6_forest_subgroups.png", 9.2, 5.4)

# ==============================================================================
# 18. Figure 7：各代谢组分年度患病率趋势
# ==============================================================================
say_sep("14. 描述性图表")

comp_long <- dat %>% filter(!is.na(met_score)) %>%
  select(year, all_of(abn_cols)) %>%
  pivot_longer(cols = all_of(abn_cols), names_to = "component", values_to = "abn") %>%
  filter(!is.na(abn)) %>%
  group_by(year, component) %>%
  summarise(prev = 100*mean(abn), n = n(), .groups = "drop") %>%
  mutate(component = recode(component,
    abn_htn="Hypertension", abn_glu="Hyperglycaemia", abn_tg="Hypertriglyceridaemia",
    abn_hdl="Low HDL-C", abn_ob="Obesity (BMI>=28)"))
p_comp <- ggplot(comp_long, aes(x = year, y = prev, colour = component, shape = component)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.6) +
  scale_y_continuous(labels = percent_format(scale = 1, suffix = "%", accuracy = 1)) +
  scale_x_continuous(breaks = 2016:2020) +
  labs(x = "Examination year", y = "Prevalence among examined participants",
       colour = NULL, shape = NULL,
       title = "Figure 7. Annual prevalence of individual metabolic abnormalities",
       subtitle = "Definitions per Chinese guidelines: SBP>=140 or DBP>=90 or treated hypertension; FBG>=6.1 mmol/L or diabetes; TG>=1.7 mmol/L; HDL-C <1.0 (men) / <1.3 (women) mmol/L; BMI>=28 kg/m2") +
  theme_pub(10.5) + theme(plot.subtitle = element_text(size = 8.2),
                          legend.text = element_text(size = 9))
save_fig(p_comp, "fig7_component_prevalence.png", 8.8, 5.4)

# ==============================================================================
# 19. Figure 8：按累积负担分组的 eGFR 纵向轨迹
# ==============================================================================
traj <- dat %>% filter(id %in% dm$id, !is.na(egfr), !is.na(met_score)) %>%
  select(id, year, egfr) %>%
  left_join(dm %>% select(id, burden_q), by = "id") %>%
  filter(!is.na(burden_q)) %>%
  group_by(year, burden_q) %>%
  summarise(mean_egfr = mean(egfr), n = n(),
            se = sd(egfr)/sqrt(n), .groups = "drop")
p_traj <- ggplot(traj, aes(x = year, y = mean_egfr, colour = burden_q, group = burden_q)) +
  geom_ribbon(aes(ymin = mean_egfr - 1.96*se, ymax = mean_egfr + 1.96*se),
              alpha = 0.16, colour = NA) +
  geom_line(linewidth = 0.95) + geom_point(size = 2.4) +
  scale_colour_manual(values = c("Q1"="#5aa469","Q2"="#8fb8de","Q3"="#e8a33d","Q4"="#d1495b"),
                      labels = c("Q1 (lowest burden)","Q2","Q3","Q4 (highest burden)")) +
  scale_x_continuous(breaks = 2016:2020) +
  labs(x = "Examination year", y = expression("Mean eGFR (mL/min/1.73 m"^2*")"),
       colour = "Cumulative burden quartile",
       title = "Figure 8. Longitudinal eGFR trajectories by cumulative burden quartile",
       subtitle = "Shaded bands = 95% CI around group means. Near-parallel trajectories are consistent with the null association reported in the primary analysis") +
  theme_pub(10.5) + theme(plot.subtitle = element_text(size = 8.6))
save_fig(p_traj, "fig8_egfr_trajectories.png", 8.4, 5.4)

# ==============================================================================
# 20. Figure 9：eGFR 斜率分布（按负担分组）
# ==============================================================================
p_slope <- ggplot(dm, aes(x = slope, fill = factor(rapid_decline))) +
  geom_histogram(binwidth = 1, colour = "white", linewidth = 0.15, alpha = 0.85) +
  geom_vline(xintercept = c(-5, -3), linetype = "dashed", colour = "grey25") +
  annotate("text", x = -5, y = Inf, label = " -5 (primary)", vjust = 1.6, hjust = -0.1, size = 3) +
  annotate("text", x = -3, y = Inf, label = " -3 (sensitivity)", vjust = 3.2, hjust = -0.1, size = 3) +
  scale_fill_manual(values = c("0"="#8fb8de","1"="#d1495b"),
                    labels = c("Stable","Rapid decline"), name = NULL) +
  coord_cartesian(xlim = quantile(dm$slope, c(0.001, 0.999), na.rm = TRUE)) +
  labs(x = expression("Individual annual eGFR slope (mL/min/1.73 m"^2*" per year)"),
       y = "Number of participants",
       title = "Figure 9. Distribution of individual eGFR slopes",
       subtitle = sprintf("Vertical dashed lines mark the rapid-decline thresholds; %.1f%% of participants fell below -5",
                          100*mean(dm$rapid_decline))) +
  theme_pub(10.5)
save_fig(p_slope, "fig9_slope_distribution.png", 8.0, 4.8)

# ==============================================================================
# 21. Figure 10：Model 2 各协变量效应森林图
# ==============================================================================
m2_plot <- r2 %>% filter(term != "(Intercept)") %>%
  mutate(label = factor(label, levels = rev(label)))
p_m2 <- ggplot(m2_plot, aes(x = OR, y = label)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = lower, xmax = upper), height = 0.2, colour = "#4a7ba7",
                linewidth = 0.7, orientation = "y") +
  geom_point(aes(size = wald_chi2), colour = "#2c6fbb") +
  geom_text(aes(x = upper, label = sprintf("%.3f (%.3f-%.3f)", OR, lower, upper)),
            hjust = -0.12, size = 2.9, colour = "grey20") +
  scale_size_continuous(range = c(2.2, 5.5), name = expression(Wald~chi^2)) +
  scale_x_log10() +
  labs(x = "Odds ratio (95% CI, log scale)", y = NULL,
       title = "Figure 10. Multivariable predictors of rapid kidney function decline (Model 2)",
       subtitle = "Point size proportional to Wald chi-square. Note the counter-intuitive positive association of baseline eGFR and the inverse association of baseline SBP") +
  theme_pub(10.5) + theme(plot.subtitle = element_text(size = 8.4))
save_fig(p_m2, "fig10_forest_model2.png", 9.0, 4.8)

# ==============================================================================
# 22. 汇总表导出
# ==============================================================================
say_sep("15. 结果汇总")

key_res <- data.frame(
  Item = c("Final analytic sample (n)", "Men, n (%)", "Median follow-up (years)",
           "Person-years of observation", "Rapid decline events, n (%)",
           "Cumulative burden, mean (SD)", "Baseline eGFR, mean (SD)",
           "Baseline SBP, mean (SD)", "Baseline BMI, mean (SD)", "Baseline age, mean (SD)",
           "Model 1 crude OR per 1-point burden",
           "Model 2 adjusted OR per 1-point burden",
           "Model 3 fully adjusted OR per 1-point burden",
           "P for trend across quartiles",
           "AUC burden alone", "AUC core model", "AUC Model 2",
           "DeLong P (core vs. burden-added)", "Continuous NRI", "IDI",
           "P_interaction sex", "P_interaction age"),
  Value = c(
    format(nrow(dm), big.mark=","),
    sprintf("%d (%.1f%%)", sum(dm$sex=="Male"), 100*mean(dm$sex=="Male")),
    as.character(median(dm$follow_years)),
    format(sum(dm$follow_years), big.mark=","),
    sprintf("%d (%.1f%%)", sum(dm$rapid_decline), 100*mean(dm$rapid_decline)),
    sprintf("%.2f (%.2f)", mean(dm$total_burden), sd(dm$total_burden)),
    sprintf("%.2f (%.2f)", mean(dm$base_egfr), sd(dm$base_egfr)),
    sprintf("%.2f (%.2f)", mean(dm$base_sbp), sd(dm$base_sbp)),
    sprintf("%.2f (%.2f)", mean(dm$base_bmi), sd(dm$base_bmi)),
    sprintf("%.2f (%.2f)", mean(dm$age), sd(dm$age)),
    fmt_or(r1, "total_burden"), fmt_or(r2, "total_burden"), fmt_or(r3, "total_burden"),
    format.pval(p_trend, eps=0.001, digits=3),
    ci_str(roc_burden), ci_str(roc_core), ci_str(roc_m2),
    format.pval(dt2$p.value, eps=0.001, digits=3),
    sprintf("%.4f (%.4f to %.4f), P=%s", ni$estimate[1], ni$lower[1], ni$upper[1],
            format.pval(ni$p[1], eps=0.001, digits=3)),
    sprintf("%.4f (%.4f to %.4f), P=%s", ni$estimate[4], ni$lower[4], ni$upper[4],
            format.pval(ni$p[4], eps=0.001, digits=3)),
    format.pval(p_int_sex, eps=0.001, digits=3),
    format.pval(p_int_age, eps=0.001, digits=3)),
  stringsAsFactors = FALSE)
write.csv(key_res, "tables/key_results_summary.csv", row.names = FALSE)
say("  [表] tables/key_results_summary.csv")
print(key_res, row.names = FALSE)

# 组分患病率年度趋势表
comp_wide <- comp_long %>% select(year, component, prev) %>%
  pivot_wider(names_from = component, values_from = prev) %>%
  mutate(across(where(is.numeric), ~ round(., 2)))
write.csv(comp_wide, "tables/component_prevalence_by_year.csv", row.names = FALSE)
say("  [表] tables/component_prevalence_by_year.csv")
print(as.data.frame(comp_wide), row.names = FALSE)

say("\n分析完成。所有图表位于 figures/，表格位于 tables/，建模数据为 analysis_data_corrected.csv")
close(LOG)
