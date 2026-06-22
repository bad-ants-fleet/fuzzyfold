use crate::Base;
use crate::parameters::parameterset::LoopEntry;

const fn convert<const N: usize>(bytes: &[u8; N]) -> [Base; N] {
    let mut out = [Base::A; N];
    let mut i = 0;
    while i < N {
        out[i] = 
            match bytes[i] {
                b'A' => Base::A,
                b'C' => Base::C,
                b'G' => Base::G,
                b'U' => Base::U,
                _ => panic!("Invalid base"),
            };
        i += 1;
    }
    out
}

macro_rules! loop_entry {
    ($seq:literal, $val:expr) => {
        LoopEntry {
            seq: &convert($seq),
            val: $val,
        }
    };
}

pub static COLLECTED_TRILOOPS_EN37: &[LoopEntry] = &[
    loop_entry!(b"CAUAG", 576), /* 574 in Thulasi et al. 2010 */ 
    loop_entry!(b"CAGAG", 554), /* 548 in Thulasi et al. 2010 */ 
    loop_entry!(b"UUUUA", 535), /* 486 in Thulasi et al. 2010 */ 
    loop_entry!(b"UAAAA", 633), /* 584 in Thulasi et al. 2010 */ 
    loop_entry!(b"CAAAG", 540), /* 538 in Thulasi et al. 2010, Groebe et al. 1988 */ 
    loop_entry!(b"CUUCG", 562), /* 556 in Thulasi et al. 2010 */
    loop_entry!(b"GUUUC", 470), /* 464 in Thulasi et al. 2010 */
    loop_entry!(b"GACAC", 644), /* 638 in Thulasi et al. 2010 */
    loop_entry!(b"UCAAA", 640), /* 591 in Thulasi et al. 2010 */
    loop_entry!(b"AUAUU", 595), /* 549 in Thulasi et al. 2010 */
    loop_entry!(b"UUAUA", 595), /* 546 in Thulasi et al. 2010 */
    loop_entry!(b"CUCCG", 523), /* 517 in Thulasi et al. 2010 */
    loop_entry!(b"GAGAC", 599), /* 593 in Thulasi et al. 2010 */
    loop_entry!(b"GAUAC", 634), /* 626 in Thulasi et al. 2010, Serra et al. 1997 */
    loop_entry!(b"GCUUC", 682), /* 676 in Thulasi et al. 2010 */
    loop_entry!(b"UACAA", 619), /* 570 in Thulasi et al. 2010 */
    loop_entry!(b"CGAAG", 465), /* 459 in Thulasi et al. 2010 */
    loop_entry!(b"CUAUG", 476), /* 470 in Thulasi et al. 2010 */
    loop_entry!(b"AAAAU", 644), /* 598 in Thulasi et al. 2010 */
    loop_entry!(b"GAAAC", 547), /* 541 in Thulasi et al. 2010 */
    loop_entry!(b"CACAG", 584), /* 578 in Thulasi et al. 2010 */
    loop_entry!(b"GCAAC", 568), /* 562 in Thulasi et al. 2010 */
    loop_entry!(b"ACAUU", 633), /* 587 in Thulasi et al. 2010 */
    loop_entry!(b"CCUUG", 453), /* 447 in Thulasi et al. 2010 */
    loop_entry!(b"UAACG", 589), /* 540 in Thulasi et al. 2010 */
    loop_entry!(b"GACCC", 587), /* 581 in Thulasi et al. 2010 */
    loop_entry!(b"UACCA", 617), /* 568 in Thulasi et al. 2010 */
    loop_entry!(b"GAAAU", 510), /* 637 in Thulasi et al. 2010, Giese et al. 1998 */

    loop_entry!(b"CAACG", 680), /* Turner et al. 2004 */
    loop_entry!(b"CUUAC", 690), /* Turner et al. 2004 */
];

pub static COLLECTED_TRILOOPS_ENTH: &[LoopEntry] = &[
    loop_entry!(b"CAUAG",   640), /* Thulasi et al. 2010 */ 
    loop_entry!(b"CAGAG",   190), /* Thulasi et al. 2010 */
    loop_entry!(b"UUUUA",   260), /* Thulasi et al. 2010 */
    loop_entry!(b"UAAAA",   890), /* Thulasi et al. 2010 */
    loop_entry!(b"CAAAG",   370), /* Thulasi et al. 2010, Groebe et al. 1988 */ 
    loop_entry!(b"CUUCG",   560), /* Thulasi et al. 2010 */
    loop_entry!(b"GUUUC",  -450), /* Thulasi et al. 2010 */
    loop_entry!(b"GACAC",   980), /* Thulasi et al. 2010 */
    loop_entry!(b"UCAAA",  1090), /* Thulasi et al. 2010 */
    loop_entry!(b"AUAUU",   510), /* Thulasi et al. 2010 */
    loop_entry!(b"UUAUA",   710), /* Thulasi et al. 2010 */
    loop_entry!(b"CUCCG",   230), /* Thulasi et al. 2010 */
    loop_entry!(b"GAGAC",   190), /* Thulasi et al. 2010 */
    loop_entry!(b"GAUAC",   -80), /* Thulasi et al. 2010, Serra et al. 1997 */
    loop_entry!(b"GCUUC",  1220), /* Thulasi et al. 2010 */
    loop_entry!(b"UACAA",   810), /* Thulasi et al. 2010 */
    loop_entry!(b"CGAAG",  -170), /* Thulasi et al. 2010 */
    loop_entry!(b"CUAUG",   -90), /* Thulasi et al. 2010 */
    loop_entry!(b"AAAAU",   870), /* Thulasi et al. 2010 */
    loop_entry!(b"GAAAC",    90), /* Thulasi et al. 2010 */
    loop_entry!(b"CACAG",   350), /* Thulasi et al. 2010 */
    loop_entry!(b"GCAAC",   590), /* Thulasi et al. 2010 */
    loop_entry!(b"ACAUU",   710), /* Thulasi et al. 2010 */
    loop_entry!(b"CCUUG", -1010), /* Thulasi et al. 2010 */
    loop_entry!(b"UAACG",   590), /* Thulasi et al. 2010 */
    loop_entry!(b"GACCC",   360), /* Thulasi et al. 2010 */
    loop_entry!(b"UACCA",   790), /* Thulasi et al. 2010 */
    loop_entry!(b"GAAAU",   330), /* Thulasi et al. 2010, Giese et al. 1998 */

    loop_entry!(b"CAACG",  2370), /* Turner et al. 2004 */
    loop_entry!(b"CUUAC",  1080), /* Turner et al. 2004 */
];

pub static COLLECTED_TETRALOOPS_EN37: &[LoopEntry] = &[
    loop_entry!(b"CGAAAG", 285), /* Sheehy et al. 2010 */
    loop_entry!(b"CGUGAG", 297), /* Sheehy et al. 2010 */
    /* loop_entry!(b"CUUCGG", 370), Turner et al. 2004 */
    loop_entry!(b"CUUCGG", 126), /* Sheehy et al. 2010 */
    loop_entry!(b"GGAAAC", 375), /* Sheehy et al. 2010 */
    loop_entry!(b"CGCAAG", 310), /* Sheehy et al. 2010 */
    loop_entry!(b"CGAGAG", 347), /* Sheehy et al. 2010 */
    loop_entry!(b"CGUAAG", 279), /* Sheehy et al. 2010 */
    loop_entry!(b"GGCAAC", 355), /* Sheehy et al. 2010 */
    loop_entry!(b"AGAAAU", 400), /* Sheehy et al. 2010 */
    loop_entry!(b"GGUGAC", 408), /* Sheehy et al. 2010 */
    loop_entry!(b"GGGGAC", 375), /* Sheehy et al. 2010 */
    loop_entry!(b"UGAAAA", 400), /* Sheehy et al. 2010 */
    /* loop_entry!(b"CUACGG", 280), Turner et al. 2004 */
    loop_entry!(b"CUACGG", 243), /* Sheehy et al. 2010 */
    loop_entry!(b"CGCGAG", 259), /* Sheehy et al. 2010 */
    loop_entry!(b"UGAAAG", 420), /* Sheehy et al. 2010 */
    loop_entry!(b"CGGAAG", 340), /* Dale et al. 2000 */
    loop_entry!(b"AAACAU", 490), /* Sheehy et al. 2010 */
    loop_entry!(b"UGAGAG", 423), /* Sheehy et al. 2010 */
    loop_entry!(b"CGAAGG", 382), /* Sheehy et al. 2010 */
    loop_entry!(b"AGCAAU", 407), /* Sheehy et al. 2010 */
    loop_entry!(b"GGAGAC", 378), /* Sheehy et al. 2010 */
    loop_entry!(b"GGGAAC", 390), /* Sheehy et al. 2010 */
    loop_entry!(b"CUAACG", 382), /* Sheehy et al. 2010 */
    loop_entry!(b"GGUAAC", 419), /* Sheehy et al. 2010 */

    loop_entry!(b"CAACGG", 550), /* Turner et al. 2004 */
    loop_entry!(b"CCAAGG", 330), /* Turner et al. 2004 */
    loop_entry!(b"CCACGG", 370), /* Turner et al. 2004 */
    loop_entry!(b"CCCAGG", 340), /* Turner et al. 2004 */
    loop_entry!(b"CCGAGG", 350), /* Turner et al. 2004 */
    loop_entry!(b"CCGCGG", 360), /* Turner et al. 2004 */
    loop_entry!(b"CCUAGG", 370), /* Turner et al. 2004 */
    loop_entry!(b"CCUCGG", 250), /* Turner et al. 2004 */
    loop_entry!(b"CUAAGG", 360), /* Turner et al. 2004 */
    loop_entry!(b"CUCAGG", 370), /* Turner et al. 2004 */
    loop_entry!(b"CUCCGG", 270), /* Turner et al. 2004 */
    loop_entry!(b"CUGCGG", 280), /* Turner et al. 2004 */
    loop_entry!(b"CUUAGG", 350), /* Turner et al. 2004 */
    loop_entry!(b"CUUUGG", 370), /* Turner et al. 2004 */
];

pub static COLLECTED_TETRALOOPS_ENTH: &[LoopEntry] = &[
    loop_entry!(b"CGAAAG",  -710), /* Sheehy et al. 2010 */
    loop_entry!(b"CGUGAG", -1000), /* Sheehy et al. 2010 */
    /* loop_entry!(b"CUUCGG", -1530), Turner et al. 2004 */
    loop_entry!(b"CUUCGG", -2100), /* Sheehy et al. 2010 */
    loop_entry!(b"GGAAAC",  -830), /* Sheehy et al. 2010 */
    loop_entry!(b"CGCAAG",  -530), /* Sheehy et al. 2010 */
    loop_entry!(b"CGAGAG",  -560), /* Sheehy et al. 2010 */
    loop_entry!(b"CGUAAG",  -810), /* Sheehy et al. 2010 */
    loop_entry!(b"GGCAAC",  -890), /* Sheehy et al. 2010 */
    loop_entry!(b"AGAAAU",  -900), /* Sheehy et al. 2010 */
    loop_entry!(b"GGUGAC",  -640), /* Sheehy et al. 2010 */
    loop_entry!(b"GGGGAC", -2280), /* Sheehy et al. 2010 */
    loop_entry!(b"UGAAAA",  -750), /* Sheehy et al. 2010 */
    /* loop_entry!(b"CUACGG", -1070), Turner et al. 2004 */
    loop_entry!(b"CUACGG", -1140), /* Sheehy et al. 2010 */
    loop_entry!(b"CGCGAG", -1190), /* Sheehy et al. 2010 */
    loop_entry!(b"UGAAAG",  -710), /* Sheehy et al. 2010 */
    loop_entry!(b"CGGAAG",  -520), /* Dale et al. 2000 */
    loop_entry!(b"AAACAU",  -240), /* Sheehy et al. 2010 */
    loop_entry!(b"UGAGAG",  -810), /* Sheehy et al. 2010 */
    loop_entry!(b"CGAAGG",  -300), /* Sheehy et al. 2010 */
    loop_entry!(b"AGCAAU",  -720), /* Sheehy et al. 2010 */
    loop_entry!(b"GGAGAC", -1060), /* Sheehy et al. 2010 */
    loop_entry!(b"GGGAAC",  -830), /* Sheehy et al. 2010 */
    loop_entry!(b"CUAACG",  -450), /* Sheehy et al. 2010 */
    loop_entry!(b"GGUAAC",  -610), /* Sheehy et al. 2010 */

    loop_entry!(b"CAACGG",   690), /* Turner et al. 2004 */
    loop_entry!(b"CCAAGG", -1030), /* Turner et al. 2004 */
    loop_entry!(b"CCACGG",  -330), /* Turner et al. 2004 */
    loop_entry!(b"CCCAGG",  -890), /* Turner et al. 2004 */
    loop_entry!(b"CCGAGG",  -660), /* Turner et al. 2004 */
    loop_entry!(b"CCGCGG",  -750), /* Turner et al. 2004 */
    loop_entry!(b"CCUAGG",  -350), /* Turner et al. 2004 */
    loop_entry!(b"CCUCGG", -1390), /* Turner et al. 2004 */
    loop_entry!(b"CUAAGG",  -760), /* Turner et al. 2004 */
    loop_entry!(b"CUCAGG",  -660), /* Turner et al. 2004 */
    loop_entry!(b"CUCCGG", -1290), /* Turner et al. 2004 */
    loop_entry!(b"CUGCGG", -1070), /* Turner et al. 2004 */
    loop_entry!(b"CUUAGG",  -620), /* Turner et al. 2004 */
    loop_entry!(b"CUUUGG",  -680), /* Turner et al. 2004 */
];

pub static COLLECTED_HEXALOOPS_EN37: &[LoopEntry] = &[
    loop_entry!(b"ACAGUACU", 280), /* Turner et al. 2004 */
    loop_entry!(b"ACAGUGAU", 360), /* Turner et al. 2004 */
    loop_entry!(b"ACAGUGCU", 290), /* Turner et al. 2004 */
    loop_entry!(b"ACAGUGUU", 180), /* Turner et al. 2004 */
];

pub static COLLECTED_HEXALOOPS_ENTH: &[LoopEntry] = &[
    loop_entry!(b"ACAGUACU", -1680), /* Turner et al. 2004 */
    loop_entry!(b"ACAGUGAU", -1140), /* Turner et al. 2004 */
    loop_entry!(b"ACAGUGCU", -1280), /* Turner et al. 2004 */
    loop_entry!(b"ACAGUGUU", -1540), /* Turner et al. 2004 */
];

