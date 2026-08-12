# Fall 2027 college applications — handoff

**Last updated:** 2026-08-10

## Primary artifact

- `fall-2027-app-timelines.html` — single-page tracker with Mermaid Gantt chart + interactive checklist (localStorage).

Open in a browser (no build step). Checklist state key: `fall2027-checklist-v10`.

**Applicant:** US citizen · CA resident · AP · **overseas Chinese (僑生)** · Fall 2027 entry.

**Checklist:** action items only. Decision dates on Gantt chart only (gold = deadline, teal = decision).

## Tracked schools (internal keys)

| Key | Display | Notes |
|-----|---------|-------|
| `uc` | UC | CA resident |
| `hku`, `cuhk`, `thu` | HK / China intl | |
| `nthu` | NTHU 僑生 direct | OGA direct admission Aug 14–Oct 15, 2026; Business Management (IBP Group B) |
| `nthui` | NTHU IBP intl portal | Optional dual pathway for 僑生; ibp.nthu.edu.tw |
| `uec` | UECFOCS | Committee bachelor's individual + joint distribution |
| `ntutw` | NTU TW intl | Foreign-national portal — 僑生 may also apply if eligible |
| `nus`, `ntu`, `um` | Singapore / Malaysia | |

**Not tracked:** NTU Taiwan 僑生 school-direct single admission — limited English-taught options; none match applicant interest.

## UECFOCS individual application — NOT graduate-only

[US & Canada area](https://cmn-hant.overseas.ncnu.edu.tw/further-study-area/national-area/america-and-canada-area/) lists separate tracks:

| Track | Who | Apply (est. Fall 2027) |
|-------|-----|------------------------|
| **學士班（個人申請）** | High school → bachelor's | Nov 1 – Dec 15, 2026; docs by Jan 6, 2027 |
| **學士班（聯合分發）** | High school → bachelor's | Feb 1 – Feb 28, 2027 |
| **研究所** | Master's / PhD | Nov 1 – Dec 15, 2026 (separate 學制 on portal) |

Portal: choose 學士班(含僑先部) vs 碩士班 vs 博士班 at [student.overseas.ncnu.edu.tw](https://student.overseas.ncnu.edu.tw/).

### NTHU via UECFOCS bachelor's individual

Per [NTHU IBP FAQ](https://ibp.nthu.edu.tw/faq-2.html): 僑生 may use committee **個人申請** for IBP **Group A + Group B** (Business Management). Also: school direct (`nthu`) and optional intl IBP portal (`nthui`). NTHU **master's** must use UECFOCS (school direct is undergrad only).

### NTU Taiwan via UECFOCS bachelor's individual

- **116學年度 (Fall 2027):** NTU 僑生 **undergrad** uses **school direct admission** on admissions.ntu.edu.tw (Aug–Sep), not the general committee individual route for most majors.
- Committee **個人申請** historically listed only niche NTU programs (e.g. 114 cycle: Japanese, Vet Med) — not the main 僑生 pathway for NTU.
- **Tracked instead:** NTU TW intl foreign portal (`ntutw`) if applying as international.

## i僑卡 (i Compatriot Card)

| System | Purpose | Required for apps? |
|--------|---------|-------------------|
| **i僑卡** | OCAC digital member card — discounts, OCAC activity auto-fill | **No** for UECFOCS/NTHU/NTU application login. **Recommended** for OCAC scholarships, work-study, activities after enrollment. [icard.taiwan-world.net](https://icard.taiwan-world.net/) — review ~3–7 business days. **Not** 僑生 identity proof. |

NTHU 僑生 direct uses 華裔身分認定初審檢核表 + residency docs, not i僑卡.

## Taiwan pathways (primary)

### NTHU — 僑生 direct (`nthu`)

[Direct admission](https://apply.nthu.edu.tw/en/article/97-overseas-chinese) — Aug 14–Oct 15, 2026. Business Management via IBP Group B.

### UECFOCS committee (`uec`)

Individual: NTHU IBP Group B yes. Joint: NTHU Group A only (not Business Management).

## Gantt milestone colors

Mermaid milestones = rotated `rect` with class `milestone`. `colorGanttMilestones()` paints `m_*` gold, `d_*` teal.

## Editing

Bump `STORAGE_KEY` when checklist IDs change.
