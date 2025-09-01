#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <cairo/cairo.h>
#include <cairo/cairo-svg.h>
#include <math.h>
#include <vector>
#include <algorithm>

#define IMAGE_SIZE      1000.0f
#define MIN_TICK_STEP   5.0f
#define MARGIN_LEFT     50.0f
#define GRAPH_MARGIN_TOP 50.0f   
#define MARGIN_TOP       150.0f  
#define MARGIN_RIGHT    100.0f
#define MARGIN_BOTTOM   50.0f
#define TEXT_HEIGHT     100.0f

#define DNA_MATCH       (1)
#define DNA_MISMATCH    (-3)
#define DNA_GAP_EXT     (2)
#define DNA_GAP_OPEN    (3)
#define DNA_GAP_FIRST   (DNA_GAP_EXT+DNA_GAP_OPEN)
#define DEBUG (0)

struct my_rect_t {
    double x0, y0, x1, y1;
    double i0, j0, i1, j1;
    double width()  { return x1 - x0; }
    double height() { return y1 - y0; }
    double toY(double i) { return ((i - i0) / (i1 - i0)) * height() + y0; }
    double toX(double j) { return ((j - j0) / (j1 - j0)) * width() + x0; }
    void incX(double x) { x0 += x; x1 += x; }
    void incY(double y) { y0 += y; y1 += y; }
    void validate() {
        if (x0 < 0) incX(-x0);
        if (y0 < 0) incY(-y0);
    }
};

struct PruningStats {
    float perc_left;
    float perc_right;
    float perc_total;
    float millions_cells;
};

static int dna_gap_open;
static int dna_gap_ext;
static int dna_match;
static int dna_mismatch;
char seq0[50], seq1[50];
int seq0_size;
int seq1_start = 0;
int seq1_end = 0;
char filetxt[256], filesvg[256], fileout[256], fileplot[256];
char filetxts[4][256], fileouts[4][256], fileplots[4][256];

float total_blocks = 0 ;
float total_pruned = 0;


static PruningStats drawPruningArea(cairo_t* cr, my_rect_t* rect) {
	cairo_save(cr);  
	//Sequence* seq0 = &job->seq[0];
	//Sequence* seq1 = &job->seq[1];
	
	const int seq0_len = seq0_size;
	const int seq1_len = (seq1_end - seq1_start);
	printf("seq1_start = %d\n", seq1_start);
	printf("seq0_len = %d, seq1_len = %d\n", seq0_len, seq1_len);

	PruningStats stats;
    stats.perc_left = 0;
    stats.perc_right = 0;
    stats.perc_total = 0;
    stats.millions_cells = 0;
	
	FILE* dump = fopen(filetxt, "rt");
	printf("Abertura feita com sucesso!\n");
	FILE* result = fopen(fileout, "wt");
	FILE* plot = fopen(fileplot, "wt");
	//fprintf(stderr, "dump %s %p\n", "pruning.txt", dump);
	if (dump == NULL){
		printf("dump não encontrado\n");
		return stats;
	} 
	if (result == NULL){
		printf("dump não encontrado\n");
		return stats;
	} 
	if (plot == NULL){
		printf("dump não encontrado\n");
		return stats;
	} 
	
	cairo_rectangle(cr, rect->x0, rect->y0, rect->width(), rect->height());
	cairo_clip(cr);
	
    std::vector<int> k_il;
	std::vector<int> k_jl;
	std::vector<int> k_ir;
	std::vector<int> k_jr;
	int d, s, e;
	int width, height;
    int count = 0;
    fread(&height, sizeof(int), 1, dump);
    fread(&width, sizeof(int), 1, dump);
    printf("Blocks: %d - %d \n", width, height);

	/*while (count < 20) {
		fread(&d, sizeof(int), 1, dump);
		fread(&s, sizeof(int), 1, dump);
		fread(&e, sizeof(int), 1, dump);
		printf("Blocks: %d - %d - %d \n", d, s ,e);
		k_d.push_back(d);
		k_e.push_back(e);
		k_s.push_back(s);
		count++;
	}*/
	int h, w, cel;
	int blocksizeH = seq1_len/width;
	int blocksizeV = seq0_len/height;
	int bpleft = 0;
	int bpright = 0;
	int nobp = 0;

    printf("height: %d - BsV: %d - width: %d - BsH: %d \n", height, blocksizeV, width, blocksizeH);
    //return;

	for (h=0; h<height; h++) {
		 for (w=0;w<width;w++) {
			int pos = h*width + w; 
			fseek(dump, (pos + 2) * sizeof(int), SEEK_SET);
			fread(&cel, sizeof(int), 1, dump);
			if (cel < 0) {  
			  //fputc('X', result);
			  if (h*blocksizeV > w*blocksizeH) { 
				  k_il.push_back(h*blocksizeV); 
				  k_jl.push_back(w*blocksizeH);
				  bpleft++;
			  }
			  else {
				  k_ir.push_back(h*blocksizeV);
				  k_jr.push_back(w*blocksizeH);
				  bpright++;
			  }
			  fprintf(plot, "%d %d\n",h*blocksizeV, w*blocksizeH); 
			  if (w != width-1)
				  fprintf(plot, "%d %d\n",h*blocksizeV, (w*blocksizeH + blocksizeH));
			  if (h != height-1)
				  fprintf(plot, "%d %d\n",(h*blocksizeV  + blocksizeV), w*blocksizeH);
			  if ((w != width-1) && (h != height-1))
				  fprintf(plot, "%d %d\n",(h*blocksizeV  + blocksizeV), (w*blocksizeH  + blocksizeH));
			}
			else {
			  //fputc(' ', result);
			  //fputc('-', result);
			nobp++;
			}
		 }
		 //fputc('\n', result);
	}

	//for (int k=1;k<k_il.size();k++)
	//	printf("%d - %d \n", k_il[k], k_jl[k]);

	float bptotal = (float)(bpleft+bpright);
	float total = (float) (bptotal + nobp);
	float cells = (float) (bptotal*blocksizeH/1000000);
	cells = cells * (float) blocksizeH;

	total_blocks += total;
	total_pruned += bptotal;


    stats.perc_left = (bpleft/bptotal)*100;
    stats.perc_right = (bpright/bptotal)*100;
    stats.perc_total = (100*bptotal/total);
    stats.millions_cells = cells;

	fprintf(result,"Height: %d - Width: %d \n", height, width);
	fprintf(result, "BPleft: %d - BPright: %d - BPtotal: %d \n", bpleft, bpright, (bpleft+bpright));
	fprintf(result, "No BP: %d -  Perc of BP: %f \n", nobp, (100*bptotal/total));
	fprintf(result, "Millions of cells pruned: %f \n", cells);
	fflush(result);
    fclose(result);
    //fclose(dump);
    //return;
	//int blocks = k_e[0];

	//int max_d = k_d[k_d.size()-1]-blocks;
	
	cairo_set_line_width (cr, 0.1);
	cairo_set_source_rgb(cr, 0.6, 0.6, 0.6);
	
	if (!k_il.empty()) {
		int i0 = k_il[0];
		int j0_global = k_jl[0] + rect->j0;
		cairo_move_to(cr, rect->toX(j0_global), rect->toY(i0));  // inicia no primeiro ponto de pruning real
		for (int i = 1; i < k_il.size(); i++) {
			int i0 = k_il[i];
			int j0_global = k_jl[i] + rect->j0;
			cairo_line_to(cr, rect->toX(j0_global), rect->toY(i0));
		}
	}


	cairo_line_to(cr, rect->toX(rect->j1), rect->toY(rect->i1));
	cairo_line_to(cr, rect->toX(rect->j0), rect->toY(rect->i1));
	// cairo_fill(cr);  

	cairo_move_to(cr, rect->toX(rect->j1), rect->toY(rect->i0));
	for (int l = 1; l < k_ir.size(); l++) {
		int i1 = k_ir[l];
		int j1_global = k_jr[l] + rect->j0;  
		cairo_line_to(cr, rect->toX(j1_global), rect->toY(i1));
	}

	cairo_line_to(cr, rect->toX(rect->j1), rect->toY(rect->i1));

	//cairo_fill(cr);
	fclose(dump);
	
	cairo_reset_clip(cr);
    cairo_restore(cr);  
    return stats;
}


static void drawBackground(cairo_t* cr, my_rect_t* rect) { 
	cairo_set_line_width (cr, 0.1);
	
	cairo_rectangle(cr, rect->x0, rect->y0, rect->width(), rect->height());
	cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
	cairo_fill(cr);
	
	cairo_rectangle(cr, rect->x0, rect->y0, rect->width(), rect->height());
	cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
	cairo_set_line_width(cr, 1);
	cairo_stroke(cr);
}

	
static void draw_accession_numbers(cairo_t* cr, my_rect_t* rect) {
	char text[500];
	cairo_text_extents_t te;
	
	cairo_set_font_size(cr, TEXT_HEIGHT/2);
	sprintf(text, "%s", seq1);
	cairo_text_extents (cr, text, &te);
	cairo_move_to(cr, rect->x0+(rect->width()-te.width)/2, rect->y1+TEXT_HEIGHT/2);
	cairo_show_text(cr, text);
	
	cairo_set_font_size(cr, TEXT_HEIGHT/2);
	sprintf(text, "%s", seq0);
	cairo_text_extents (cr, text, &te);
	cairo_move_to(cr, rect->x1+TEXT_HEIGHT/6, rect->y0+(rect->height()-te.width)/2);
	cairo_save (cr);
	cairo_rotate (cr, M_PI / 2.0);
	cairo_show_text(cr, text);
	cairo_restore (cr);
}	
	
static void draw_ticks(cairo_t* cr, my_rect_t* rect) {
	int diff_i = rect->i1-rect->i0;
	int diff_j = rect->j1-rect->j0;
	
	int step_ticks_i = diff_i/3;
	int step_ticks_j = diff_j/3;
	step_ticks_i = pow(10, (int)(log10(step_ticks_i)));
	step_ticks_j = pow(10, (int)(log10(step_ticks_j)));
	if ((rect->toX(step_ticks_j)-rect->toX(0))/10 < MIN_TICK_STEP) {
		step_ticks_j *= 10;
	}
	if ((rect->toY(step_ticks_i)-rect->toY(0))/10 < MIN_TICK_STEP) {
		step_ticks_i *= 10;
	}
	int tick_i0 = (rect->i0/(step_ticks_i/10)+1)*(step_ticks_i/10);
	int tick_j0 = (rect->j0/(step_ticks_j/10)+1)*(step_ticks_j/10);
	
	for (int i=tick_i0; i<=rect->i1; i+=step_ticks_i/10) {
		int len = (i%(step_ticks_i))==0?TEXT_HEIGHT/3:TEXT_HEIGHT/8;
		cairo_move_to(cr, rect->x1-len, rect->toY(i));
		cairo_line_to(cr, rect->x1, rect->toY(i));
		cairo_stroke(cr);
		cairo_move_to(cr, rect->x0, rect->toY(i));
		cairo_line_to(cr, rect->x0+len, rect->toY(i));
		cairo_stroke(cr);
	}
	for (int j=tick_j0; j<=rect->j1; j+=step_ticks_j/10) {
		int len = (j%(step_ticks_j))==0?TEXT_HEIGHT/3:TEXT_HEIGHT/8;
		cairo_move_to(cr, rect->toX(j), rect->y1);
		cairo_line_to(cr, rect->toX(j), rect->y1-len);
		cairo_stroke(cr);
		cairo_move_to(cr, rect->toX(j), rect->y0);
		cairo_line_to(cr, rect->toX(j), rect->y0+len);
		cairo_stroke(cr);
	}
	
	char text[50];
	cairo_text_extents_t te;

	cairo_set_font_size(cr, TEXT_HEIGHT/3);
	sprintf(text, "%.0f", rect->j0);
	cairo_move_to(cr, rect->x0, rect->y0-TEXT_HEIGHT/12);
	cairo_show_text(cr, text);
	sprintf(text, "%.0f", rect->j1);
	cairo_text_extents (cr, text, &te);
	cairo_move_to(cr, rect->x1-te.width, rect->y0-TEXT_HEIGHT/12);
	cairo_show_text(cr, text);
	
	
	cairo_set_font_size(cr, TEXT_HEIGHT/3);
	sprintf(text, "%.0f", rect->i0);
	cairo_move_to(cr, rect->x0-TEXT_HEIGHT/3, rect->y0);
	cairo_save (cr);
	cairo_rotate (cr, M_PI / 2.0);
	cairo_show_text(cr, text);
	cairo_restore(cr);
	
	sprintf(text, "%.0f", rect->i1);
	cairo_text_extents (cr, text, &te);
	cairo_move_to(cr, rect->x0-TEXT_HEIGHT/3, rect->y1-te.width);
	cairo_save (cr);
	cairo_rotate (cr, M_PI / 2.0);
	cairo_show_text(cr, text);
	cairo_restore(cr);
	
}

static void drawPruningStats(cairo_t* cr, my_rect_t* rect, const PruningStats& stats, int index) {
    char text[256];
    cairo_text_extents_t te;

	float pruning_total = 0;
    
    cairo_set_source_rgb(cr, 1.0, 0.0, 0.0); 
    cairo_set_font_size(cr, TEXT_HEIGHT/3.5);
    
    double x = rect->x1 -285;
    double y = rect->y0 -100;
    
	switch (index) {
		case 0:
			sprintf(text, "GPU 1_Left Stats:");
			break;
		case 1:
			sprintf(text, "GPU 1_Right Stats:");
			break;
		case 2:
			sprintf(text, "GPU 2_Left Stats:");
			break;
		case 3:
			sprintf(text, "GPU 2_Right Stats:");
			break;
		default:
			return;
	}

	// Desenha título
	cairo_move_to(cr, x, y);
	cairo_show_text(cr, text);

	// Se for GPU 2_Right, também desenha soma
	if (index == 3) {
		pruning_total = (total_pruned / total_blocks) * 100;
		cairo_move_to(cr, x + 200, y+50); // desloca horizontalmente
		sprintf(text, "Sum: %.1f%%", pruning_total);
		cairo_show_text(cr, text);
	}

	y += 40;

    /*sprintf(text, "     Left: %.1f%%", stats.perc_left);
    cairo_move_to(cr, x, y);
    cairo_show_text(cr, text);
    y += 40; */
    
    /*sprintf(text, "     Right: %.1f%%", stats.perc_right);
    cairo_move_to(cr, x, y);
    cairo_show_text(cr, text);
    y += 40;*/
    
    sprintf(text, "Total: %.1f%%", stats.perc_total);
    cairo_move_to(cr, x, y);
    cairo_show_text(cr, text);

    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0); 

}

void drawAlignmentVerticalSlices(int widths[4]) {
    cairo_surface_t *surface;
    cairo_t *cr;

    int total_width = widths[0] + widths[1] + widths[2] + widths[3];
    float size_y = IMAGE_SIZE;
    float size_x = IMAGE_SIZE * ((float)total_width / seq0_size);

    surface = cairo_svg_surface_create(filesvg, size_x + MARGIN_LEFT + MARGIN_RIGHT, size_y + MARGIN_TOP + MARGIN_BOTTOM);
    cr = cairo_create(surface);

    cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
    cairo_rectangle(cr, 0.0, 0.0, size_x + MARGIN_LEFT + MARGIN_RIGHT, size_y + MARGIN_TOP + MARGIN_BOTTOM);
    cairo_fill(cr);
	printf("size_x = %.2f, size_y = %.2f\n", size_x, size_y); 
    int j0 = 0;
    float px0 = MARGIN_LEFT;
    for (int k = 0; k < 4; ++k) {
        int j1 = j0 + widths[k];
        float px1 = px0 + ((float)widths[k] / total_width) * size_x;

        my_rect_t rect;
        rect.i0 = 0; rect.i1 = seq0_size;
        rect.j0 = j0; rect.j1 = j1;
        rect.x0 = px0; rect.x1 = px1;
        rect.y0 = MARGIN_TOP; rect.y1 = MARGIN_TOP + size_y;

        strcpy(filetxt, filetxts[k]);
        strcpy(fileout, fileouts[k]);
        strcpy(fileplot, fileplots[k]);

        seq1_start = j0;
        seq1_end = j1;
		printf("Faixa %d: j0=%d, j1=%d, px0=%.2f, px1=%.2f\n", k, j0, j1, px0, px1);

	    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0); 
		cairo_set_line_width(cr, 2.0);
		cairo_rectangle(cr, rect.x0, rect.y0, rect.width(), rect.height());
		cairo_stroke(cr);
		
        PruningStats stats = drawPruningArea(cr, &rect);
        drawPruningStats(cr, &rect, stats, k);

        j0 = j1;
        px0 = px1;
    }

    my_rect_t full_rect;
    full_rect.i0 = 0; full_rect.i1 = seq0_size;
    full_rect.j0 = 0; full_rect.j1 = j0;
    full_rect.x0 = MARGIN_LEFT; full_rect.x1 = px0;
    full_rect.y0 = MARGIN_TOP; full_rect.y1 = MARGIN_TOP + size_y;
    draw_ticks(cr, &full_rect);

    cairo_destroy(cr);
    cairo_surface_destroy(surface);
}


int main(int argc, char** argv) {
    if (argc < 10) {
        fprintf(stderr, "Uso: %s f1 f2 f3 f4 w1 w2 w3 w4 seq0_size\n", argv[0]);
        return 1;
    }

    dna_gap_open = DNA_GAP_OPEN;
    dna_gap_ext  = DNA_GAP_EXT;
    dna_match    = DNA_MATCH;
    dna_mismatch = DNA_MISMATCH;

    const char* files_in[4] = { argv[1], argv[2], argv[3], argv[4] };
    int widths[4] = { atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8]) };
    seq0_size = atoi(argv[9]);

    for (int i = 0; i < 4; ++i) {
        sprintf(filetxts[i], "%s.txt", files_in[i]);
        sprintf(fileouts[i], "%s.out", files_in[i]);
        sprintf(fileplots[i], "%s.prn", files_in[i]);
    }

    sprintf(filesvg, "%s.svg", files_in[0]);

    drawAlignmentVerticalSlices(widths);
    return 0;
}
