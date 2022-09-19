#!/usr/bin/env scheme-script

(import (cslib all))

(define path-data
  (make-parameter #f))

(define all-data
  (make-parameter #f))

(define total-site-list
  (make-parameter #f))

(define lam-list
  (make-parameter #f))

(define (get-all-data)
  (map (lambda (path-total-site)
         (cons
           (list 'total-site (scpair-find "total_site=" (make-scpair path-total-site)))
           (map (lambda (path-lam)
                  (cons
                    (list 'lam (scpair-find-f string->number "lambda=" (make-scpair path-lam)))
                    (map load-obj (glob-expand-at path-lam"mass_sqr=*"))))
                (glob-expand-at path-total-site "lambda=*"))))
       (glob-expand-at (path-data) "total_site=*")))

(define (get-total-site-list)
  (map (lambda (d) (car (rec-lookup d 'total-site))) (all-data)))

(define (get-lam-list)
  (list-sort < (map (lambda (d) (car (rec-lookup d 'lam))) (cdr (car (all-data))))))

(define (get-mass-sqr-and-m-eff v)
  ; mass-sqr m-eff-sqr m-eff-sqr-err m-eff m-eff-err
  (let* ([mass-sqr (car (rec-lookup v 'mass-sqr))]
         [m-eff-list (rec-lookup v 'n-block-32 'm_eff)]
         [m-eff (list-ref m-eff-list 0)]
         [m-eff-err (list-ref m-eff-list 1)])
    (list mass-sqr (sqr m-eff) (* 2 m-eff m-eff-err) m-eff m-eff-err)))

(define (get-mass-sqr-and-m-eff-datatable vs)
  (table->datatable
    (list-sort
      (on < car)
      (map get-mass-sqr-and-m-eff vs))))

(define (get-mass-sqr-and-m-eff-and-v-eff v)
  ; mass-sqr m-eff m-eff-err v-eff v-eff-err
  (let* ([mass-sqr (car (rec-lookup v 'mass-sqr))]
         [m-eff-list (rec-lookup v 'n-block-32 'm_eff)]
         [v-eff-list (rec-lookup v 'n-block-32 'v_eff)])
    (append (list mass-sqr) m-eff-list v-eff-list)))

(define (get-mass-sqr-and-m-eff-and-v-eff-datatable vs)
  (table->datatable
    (list-sort
      (on < car)
      (map get-mass-sqr-and-m-eff-and-v-eff vs))))

(define (mk-plot-mass-dep total-site lam)
  (let* ([vs (rec-lookup (all-data) (list 'total-site total-site) (list 'lam lam))]
         [dt (get-mass-sqr-and-m-eff-datatable vs)])
    (plot-save
      (format "plots/lambda=~a/total-site=~a/mass-dep.eps.pdf.png" lam total-site)
      (cons "table.txt" dt)
      "set size 2.0,2.0"
      "set xtics 0.1"
      (format "set title '$\\lambda=~a$'" lam)
      "set xlabel '$m^2$'"
      "set ylabel '$m_\\text{eff}^2$'"
      (mk-plot-line
        "plot [:] [:]"
        "0 not"
        (format
          "'table.txt' u 1:2:3 w yerrorb t '~a'" total-site)))))

(define (mk-plot-potential-dep total-site lam)
  (let* ([vs (rec-lookup (all-data) (list 'total-site total-site) (list 'lam lam))]
         [dt (get-mass-sqr-and-m-eff-and-v-eff-datatable vs)])
    (plot-save
      (format "plots/lambda=~a/total-site=~a/potential-dep.eps.pdf.png" lam total-site)
      (cons "table.txt" dt)
      "set size 2.0,2.0"
      "set xtics 0.1"
      (format "set title '$\\lambda=~a$'" lam)
      "set xlabel '$m_\\text{eff}$'"
      "set ylabel '$V_\\text{eff}$'"
      (mk-plot-line
        "plot [:] [:]"
        "0 not"
        (format "'table.txt' u 2:4:3:5 w xyerrorb t '~a'" total-site)))))

(define (mk-plot-potential-dep-cmp-lam total-site lam1 lam2)
  (let* ([vs1 (rec-lookup (all-data) (list 'total-site total-site) (list 'lam lam1))]
         [vs2 (rec-lookup (all-data) (list 'total-site total-site) (list 'lam lam2))]
         [dt1 (get-mass-sqr-and-m-eff-and-v-eff-datatable vs1)]
         [dt2 (get-mass-sqr-and-m-eff-and-v-eff-datatable vs2)])
    (plot-save
      (format "plots/lambda=~a/total-site=~a/potential-dep-cmp-lam.eps.pdf.png" lam1 total-site)
      (cons "table-1.txt" dt1)
      (cons "table-2.txt" dt2)
      "set size 2.0,2.0"
      "set xtics 0.1"
      (format "set title '~a'" total-site)
      "set xlabel '$m_\\text{eff}$'"
      "set ylabel '$V_\\text{eff}$'"
      (mk-plot-line
        "plot [:] [:]"
        "0 not"
        (format "'table-1.txt' u 2:4:3:5 w xyerrorb t '$\\lambda=~a$'" lam1)
        (format "'table-2.txt' u 2:4:3:5 w xyerrorb t '$\\lambda=~a$'" lam2)))))

(define (mk-plot-potential-dep-cmp-4-8 lam)
  (let* ([total-site-4 "4x4x4x256"]
         [total-site-8 "8x8x8x512"]
         [vs-4 (rec-lookup (all-data) (list 'total-site total-site-4) (list 'lam lam))]
         [vs-8 (rec-lookup (all-data) (list 'total-site total-site-8) (list 'lam lam))]
         [dt-4 (get-mass-sqr-and-m-eff-and-v-eff-datatable vs-4)]
         [dt-8 (get-mass-sqr-and-m-eff-and-v-eff-datatable vs-8)])
    (plot-save
      (format "plots/lambda=~a/cmp-4-8/potential-dep.eps.pdf.png" lam)
      (cons "table-4.txt" dt-4)
      (cons "table-8.txt" dt-8)
      "set size 2.0,2.0"
      "set xtics 0.1"
      (format "set title '$\\lambda=~a$'" lam)
      "set xlabel '$m_\\text{eff} L$'"
      "set ylabel '$V_\\text{eff} L$'"
      (mk-plot-line
        "plot [0.2:1.6] [:]"
        "0 not"
        (format "'table-4.txt' u (4*$2):(4*$4):(4*$3):(4*$5) w xyerrorb t '~a'" total-site-4)
        (format "'table-8.txt' u (8*$2):(8*$4):(8*$3):(8*$5) w xyerrorb t '~a'" total-site-8)))))

(define (mk-plot-mass-dep-all)
  (for-each
    (lambda (total-site)
      (for-each
        (lambda (lam)
          (fork-exec
            (mk-plot-mass-dep total-site lam)))
        (lam-list)))
    (total-site-list)))

(define (mk-plot-potential-dep-all)
  (for-each
    (lambda (total-site)
      (for-each
        (lambda (lam)
          (fork-exec
            (mk-plot-potential-dep total-site lam)))
        (lam-list)))
    (total-site-list)))

(define (mk-plot-potential-dep-cmp-lam-all)
  (for-each
    (lambda (total-site)
      (for-each
        (lambda (lam1 lam2)
          (fork-exec
            (mk-plot-potential-dep-cmp-lam total-site lam1 lam2)))
        (init (lam-list)) (tail (lam-list))))
    (total-site-list)))

(define (mk-plot-potential-dep-cmp-4-8-all)
  (for-each
    (lambda (lam)
      (fork-exec
        (mk-plot-potential-dep-cmp-4-8 lam)))
    (lam-list)))

(print "hello")

(fork-limit 1)

(path-data "results")

(all-data (get-all-data))

(total-site-list (get-total-site-list))

(lam-list (get-lam-list))

(print (total-site-list) (lam-list))

(print (get-mass-sqr-and-m-eff-datatable (rec-lookup (all-data) (list 'total-site "4x4x4x256") (list 'lam 32.0)))) 

(print (get-mass-sqr-and-m-eff-and-v-eff-datatable (rec-lookup (all-data) (list 'total-site "4x4x4x256") (list 'lam 32.0)))) 

(mk-plot-potential-dep-cmp-4-8-all)

(mk-plot-potential-dep-cmp-lam-all)

(mk-plot-potential-dep-all)

(mk-plot-mass-dep-all)

(wait-all)
